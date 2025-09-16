#include <biovoltron/algo/align/tailor/tailor.hpp>
#include <biovoltron/algo/annotate/annotator.hpp>
#include <biovoltron/utility/threadpool/threadpool.hpp>
#include <biovoltron/utility/istring.hpp>
#include <biovoltron/file_io/fasta.hpp>
#include <biovoltron/file_io/fastq.hpp>
#include <biovoltron/file_io/gff.hpp>
#include <catch.hpp>

#include <experimental/random>
#include <ranges>
#include <filesystem>
#include <fstream>
#include <iostream>//debug
#include <chrono>
#include <mutex>
#include <shared_mutex>
#include <list>

namespace bio = biovoltron;
namespace ranges = std::ranges;
const auto data_path = std::filesystem::path{DATA_PATH};

auto get_rand_seq(std::uint32_t base_per_chrom) {
  auto gen_rand_base = []() { 
    return bio::Codec::to_char(std::experimental::randint(0, 3)); 
  };
  auto seq = std::string(base_per_chrom, 'A');
  ranges::generate(seq, gen_rand_base);
  return seq;
}

struct IOQueue: std::shared_mutex{
  void push_job(bio::Alignment& aln){
    std::unique_lock lock(*this);
    aln_queue.emplace_back(aln);
  }
  void load_job(std::list<bio::Alignment>& dest){
    std::unique_lock lock(*this);
    dest.splice(dest.end(), aln_queue);
  }
  inline bool all_completed(){
    std::shared_lock lock(*this);
    return no_job_will_come and aln_queue.empty();
  }
  inline bool queue_empty(){
    std::shared_lock lock(*this);
    return aln_queue.empty();
  }
  void finish_thread(){
    std::unique_lock lock(*this);
    this->no_job_will_come = true;
  }
private:
  bool no_job_will_come = false;
  std::list<bio::Alignment> aln_queue;
};

TEST_CASE("Generate tailor pipeline data") {
  std::experimental::reseed(6); // Fix rand seed for easy debugging.

  constexpr auto CHROM_NUM = 3;
  constexpr auto BASE_PER_CHROM = 1000;
  constexpr auto READ_LEN = 25;
  constexpr auto FEATURE_LEN = 60;

  // Generate reference.
  auto ref = std::vector<bio::FastaRecord<false>>{};
  for (auto i = 1; i < CHROM_NUM+1; i++)
    ref.emplace_back("chr" + std::to_string(i), get_rand_seq(BASE_PER_CHROM));

  // Generate multimap reads.
  auto multimap_reads = std::vector<bio::FastqRecord<>>(20);
  for (auto i = 0, j = 0; i < multimap_reads.size(); i++) {
    auto& chrom = ref[i % CHROM_NUM];
    auto seq = chrom.seq.substr(j, READ_LEN);
    chrom.seq.replace(j + READ_LEN + 5, READ_LEN, seq);
    multimap_reads[i].seq = seq;
    multimap_reads[i].qual = std::string(READ_LEN, '!');
    if (i < multimap_reads.size() / 2)
      multimap_reads[i].name = "a-multi" + std::to_string(i);
    else 
      multimap_reads[i].name = "b-multi" + std::to_string(i);
    if (i % CHROM_NUM == CHROM_NUM - 1) j += 3*READ_LEN;
  }
  // Check every multimap read have at least 2 mapping.
  for (auto i = 0; i < multimap_reads.size(); i++) {
    const auto& chrom = ref[i%CHROM_NUM];
    const auto pos = chrom.seq.find(multimap_reads[i].seq);
    REQUIRE(pos != std::string::npos);
    REQUIRE(chrom.seq.find(multimap_reads[i].seq, pos + multimap_reads[i].seq.size()) != std::string::npos);
  }

  // Generate reverse complement reference.
  auto rc_ref = std::vector<bio::FastaRecord<false>>{};
  ranges::transform(ref, std::back_inserter(rc_ref), 
    [](auto record) { 
      record.seq = bio::Codec::rev_comp(record.seq);
      return record;
  });

  // Generate unmap reads.
  auto unmap_reads = std::vector<bio::FastqRecord<>>(20);
  auto contains = [](const std::string& read) {
    return [read](const std::string& seq) { return seq.find(read) != std::string::npos; };
  };
  for (auto i = 0; i < unmap_reads.size(); i++) {
    do unmap_reads[i].seq = get_rand_seq(READ_LEN);
    while (ranges::any_of(ref, contains(unmap_reads[i].seq), &bio::FastaRecord<false>::seq));
    unmap_reads[i].qual = std::string(READ_LEN, '!');
    if (i < unmap_reads.size() / 2)
      unmap_reads[i].name = "a-un" + std::to_string(i);
    else 
      unmap_reads[i].name = "b-un" + std::to_string(i);
  }
  // Check you can't find unmap reads in ref.
  for (const auto& record : unmap_reads)
    for (const auto& chrom : ref) 
      REQUIRE(chrom.seq.find(record.seq) == std::string::npos);

  // Generate features.
  auto feats = std::vector<bio::GffRecord>(10, bio::GffRecord{.source="Human", .type="gene"});
  for (auto i = 0; i < feats.size(); i++) {
    feats[i].attrs = "ID=gene" + std::to_string(i);

    if (i % 2 == 0) feats[i].strand = '+';
    else feats[i].strand = '-';

    feats[i].seqid = std::to_string((i%CHROM_NUM) + 1);
    feats[i].start = std::experimental::randint(0, BASE_PER_CHROM - FEATURE_LEN);
    feats[i].end = feats[i].start + FEATURE_LEN - 1;
  }
  // Check features unique.
  {
    auto unique_feats = decltype(feats){};
    ranges::unique_copy(feats, std::back_inserter(unique_feats));
    REQUIRE(feats.size() == unique_feats.size());
  }

  // Generate unique reads.
  auto unique_reads = std::vector<bio::FastqRecord<>>(160);
  for (auto i = 0; i < unique_reads.size(); i++) {
    unique_reads[i].qual = std::string(READ_LEN, '!');
    if (i < unique_reads.size() / 2)
      unique_reads[i].name = "a-unique" + std::to_string(i);
    else 
      unique_reads[i].name = "b-unique" + std::to_string(i);
  }
  for (auto idx = 0; idx < unique_reads.size();) {
    for (const auto& feat : feats | std::views::take(5)) {
      const auto num = (feat.strand == '+') ? 20 : 10;
      const auto iv = bio::Interval{feat};
      const auto chrom_idx = std::stoi(iv.chrom) - 1;
      const auto& seq = (feat.strand == '+') ? ref[chrom_idx].seq : rc_ref[chrom_idx].seq;
      for (auto i = 0; i < num; i++, idx++)
        unique_reads[idx].seq = seq.substr(iv.begin + i, READ_LEN);
    }
  }
  // Check ref and rc_ref contains one and only one uniqe reads.
  for (const auto& read : unique_reads) {
    auto cnt = 0;
    for (const auto& chrom : ref) {
      auto pos = chrom.seq.find(read.seq);
      if (pos != std::string::npos) {
        cnt++;
        REQUIRE(chrom.seq.find(read.seq, pos + READ_LEN) == std::string::npos);
      }
    }
    for (const auto& chrom : rc_ref) {
      auto pos = chrom.seq.find(read.seq);
      if (pos != std::string::npos) {
        cnt++;
        REQUIRE(chrom.seq.find(read.seq, pos + READ_LEN) == std::string::npos);
      }
    }
    REQUIRE(cnt == 1);
  }

  // Build index.
  {
    auto index = bio::Index{5};
    index.make_index(ref);
    auto ofs = std::ofstream{"ref.idx"};
    index.save(ofs);
  }
  {
    auto index = bio::Index{5};
    index.make_index(rc_ref);
    auto ofs = std::ofstream{"rc_ref.idx"};
    index.save(ofs);
  }
  {
    auto ofs = std::ofstream{"ref.fa"};
    for (const auto& record : ref)
      ofs << record << '\n';
  }
  {
    auto ofs = std::ofstream{"rc_ref.fa"};
    for (const auto& record : rc_ref)
      ofs << record << '\n';
  }
  {
    auto ofs = std::ofstream{"ref.gff"};
    for (auto feat : feats){
      if (feat.strand == '-'){
        // gff define that start & end pos should be on '+' strand no matter feat on which strand
        int tmp = BASE_PER_CHROM - feat.start - 1;
        feat.start = BASE_PER_CHROM - feat.end - 1;
        feat.end = tmp;
      }
      ofs << feat << '\n';
    }
  }
  {
    auto ofs = std::ofstream{"a.fq"};
    for (const auto& record : unmap_reads | std::views::take(static_cast<int>(unmap_reads.size() / 2)))
      ofs << record << '\n';
    for (const auto& record : multimap_reads | std::views::take(static_cast<int>(multimap_reads.size() / 2)))
      ofs << record << '\n';
    for (const auto& record : unique_reads | std::views::take(static_cast<int>(unique_reads.size() / 2)))
      ofs << record << '\n';
  }
  {
    auto ofs = std::ofstream{"b.fq"};
    for (const auto& record : unmap_reads | std::views::drop(static_cast<int>(unmap_reads.size() / 2)))
      ofs << record << '\n';
    for (const auto& record : multimap_reads | std::views::drop(static_cast<int>(multimap_reads.size() / 2)))
      ofs << record << '\n';
    for (const auto& record : unique_reads | std::views::drop(static_cast<int>(unique_reads.size() / 2)))
      ofs << record << '\n';
  }
}

TEST_CASE("Generate parallel tailor data") {
  std::experimental::reseed(42); // Fix rand seed for easy debugging.

  constexpr auto CHROM_NUM = 1;
  constexpr auto BASE_PER_CHROM = 10000000;
  constexpr auto READ_LEN = 19; // should greater than tailor seed_len
  constexpr auto FEATURE_LEN = 50;
  constexpr auto READ_COUNT = 1000000;
  constexpr auto FEATURE_COUNT = 10;
  
  // Generate reference.
  auto ref = std::vector<bio::FastaRecord<false>>{};
  for (auto i = 1; i < CHROM_NUM+1; i++)
    ref.emplace_back("chr" + std::to_string(i), get_rand_seq(BASE_PER_CHROM));
  
  // Generate reverse complement reference.
  auto rc_ref = std::vector<bio::FastaRecord<false>>{};
  ranges::transform(ref, std::back_inserter(rc_ref), 
    [](auto record) { 
      record.seq = bio::Codec::rev_comp(record.seq);
      return record;
  });
  // Generate index
  auto index = bio::Index{5};
  index.make_index(ref);
  auto rc_index = bio::Index{5};
  rc_index.make_index(rc_ref);
  
  // string::find is too slow for large ref
  auto search = [&index = std::as_const(index)](const std::string& str_query){
    auto istr_query = bio::Codec::to_istring(str_query);
    auto query = bio::istring_view(istr_query);
    auto beg = size_t{};
    auto end = index.get_bwt_size();

    while (!query.empty()) {
      if (end <= beg) break;
      beg = index.lf_mapping(query.back(), beg);
      end = index.lf_mapping(query.back(), end);
      query.remove_suffix(1);
    }
    return end - beg;
  };

  // Generate features.
  auto feats = std::vector<bio::GffRecord>(FEATURE_COUNT, bio::GffRecord{.source="Human", .type="gene"});
  for (auto i = 0; i < feats.size(); i++) {
    feats[i].attrs = "ID=gene" + std::to_string(i);
    feats[i].strand = '+';
    feats[i].seqid = "chr1";
    feats[i].start = std::experimental::randint(0, BASE_PER_CHROM - FEATURE_LEN);
    feats[i].end = feats[i].start + FEATURE_LEN - 1;
  }
  // Check features unique.
  {
    auto unique_feats = decltype(feats){};
    ranges::unique_copy(feats, std::back_inserter(unique_feats));
    REQUIRE(feats.size() == unique_feats.size());
  }
  auto all_feat_names = std::vector<std::string>{};
  auto genes = bio::Annotator<bio::GffRecord>{};
  {
    auto gff_ofs = std::ofstream{"parallel.gff"};
    for(const auto& feat: feats){
      genes.insert(feat);
      all_feat_names.emplace_back(feat.attrs.substr(3));
      gff_ofs << feat << '\n';
    }
  }
  genes.index();

  std::vector<std::string> samples = {std::string{"parallel_1.fq"}, std::string{"parallel_2.fq"}};
  auto answer = std::map<std::string, std::vector<int>>{
    {samples[0], {635705, 23, 364272}},
    {samples[1], {637032, 18, 362950}}
  };
  for(auto& sample: samples){
    auto fq_ofs = std::ofstream{sample};

    // initialize gene map
    auto gene_map_count = std::map<std::string, int>();
    for(auto& feat_name: all_feat_names)
      gene_map_count[feat_name] = 0;
    
    auto uniq_cnt = 0, multi_cnt = 0, unmap_cnt = 0;
    for(auto i = 0; i < READ_COUNT; ++i) {
      auto record = bio::FastqRecord<false>{};
      uint32_t start_pos = std::experimental::randint(0, BASE_PER_CHROM - READ_LEN); // buffer for features
      auto read = ref[0].seq.substr(start_pos, READ_LEN);
      if(std::experimental::randint(0,10) >= 7) 
        record.seq = std::string(READ_LEN, 'A') + read;
      else
        record.seq = read;
      record.qual = std::string(record.seq.length(), '!');

      auto hit_cnt = search(record.seq);
      if(hit_cnt == 1) {
        // uniq map reads
        record.name = std::string{"uniq-"} + std::to_string(i);
        auto res = genes.find(bio::Interval{"chr1", start_pos, start_pos + record.seq.length(), '+'});
        if(res.size() == 1){
          auto gene_name = res.front().attrs.substr(3);
          ++gene_map_count[gene_name];
        }
        ++uniq_cnt;
      } else if(hit_cnt > 1) {
        // multi map reads
        record.name = std::string{"multi-"} + std::to_string(i);
        ++multi_cnt;
      } else {
        // unmap reads
        record.name = std::string{"unmap-"} + std::to_string(i);
        ++unmap_cnt;
      }
      fq_ofs << record << '\n';
    }
    CHECK(uniq_cnt == answer[sample][0]);
    CHECK(multi_cnt == answer[sample][1]);
    CHECK(unmap_cnt == answer[sample][2]);
  }

  // Output index
  {
    auto ofs = std::ofstream{"parallel_ref.idx"};
    index.save(ofs);
  }
  {
    auto ofs = std::ofstream{"parallel_rc_ref.idx"};
    rc_index.save(ofs);
  }
}

TEST_CASE("Tailor pipeline") {
  // Prepair index.
  auto index = bio::Index{};
  {
    auto ifs = std::ifstream{"ref.idx"};
    index.load(ifs);
  }
  auto rc_index = bio::Index{};
  {
    auto ifs = std::ifstream{"rc_ref.idx"};
    rc_index.load(ifs);
  }

  // Init tailor.
  auto tailor = bio::Tailor{index, rc_index};
  tailor.allow_seed_mismatch = true;

  // Get raw count expression matrix. 
  auto expr_mat_b = std::unordered_map<std::string, std::uint32_t>{};
  auto expr_mat_a = std::unordered_map<std::string, std::uint32_t>{};
  // Load features into annotator.
  auto genes = bio::Annotator<bio::GffRecord>{};
  {
    auto ifs = std::ifstream{"ref.gff"};
    auto record = bio::GffRecord{};
    while (ifs >> record) {
      // Note: Add prefix "chr" to seqname.
      if (record.seqid[0] != 'c')
        record.seqid = "chr" + record.seqid;
      genes.insert(record);
    }
  }
  genes.index(); // Note: !!!!!Important remember to index.

  // Pipeline
  {
    // Sample alignment.
    auto aligner = [&tailor](const auto& record){
      return tailor.search(record);
    };

    // Output SAM file.
    auto ofs = std::ofstream{"a.sam"};
    auto sam_writer = [&ofs](auto&& aln){
      ranges::for_each(
        aln_to_sam_list(aln.first),
        [&ofs](auto&& sam_record){
          ofs << sam_record << '\n';
        }
      );
      ranges::for_each(
        aln_to_sam_list(aln.second),
        [&ofs](auto&& sam_record){
          ofs << sam_record << '\n';
        }
      );
      return std::forward<decltype(aln)>(aln);
    };

    // Matching gene interval
    auto annotator = [&genes](auto&& aln){
      bool hit_f = aln.first.hits.size() == 1;
      bool hit_r = aln.second.hits.size() == 1;

      if (hit_f && hit_r) {
        auto merged = genes.find(aln.first.hits.front().intv);
        auto result_r = genes.find(aln.second.hits.front().intv);
        merged.insert(merged.end(), result_r.begin(), result_r.end());
        return merged;
      }

      if (hit_f) {
        return genes.find(aln.first.hits.front().intv);
      }

      if (hit_r) {
        return genes.find(aln.second.hits.front().intv);
      }

      return std::vector<bio::GffRecord>{};
    };

    auto ifs = std::ifstream{"a.fq"};
    auto uniq_results = ranges::istream_view<bio::FastqRecord<>>(ifs)
      | ranges::views::transform(aligner)
      | ranges::views::transform(sam_writer)
      | ranges::views::transform(annotator);

    ranges::for_each(
      uniq_results,
      [&expr_mat_a](const auto& results){
        // Discard ambiguous annotation.
        if(results.size() == 1){
          const auto& gff_record = results.front();
          // Note: A ugly way to get gene name in attrs of gff record.
          const auto gene_name = gff_record.attrs.substr(3);
          if (expr_mat_a.contains(gene_name)) expr_mat_a[gene_name] += 1;
          else expr_mat_a[gene_name] = 1;
        }
      }
    );
    ifs = std::ifstream{"b.fq"};
    ofs = std::ofstream{"b.sam"};
    uniq_results = ranges::istream_view<bio::FastqRecord<>>(ifs)
      | ranges::views::transform(aligner)
      | ranges::views::transform(sam_writer)
      | ranges::views::transform(annotator);

    ranges::for_each(
      uniq_results,
      [&expr_mat_b](const auto& results){
        // Discard ambiguous annotation.
        if(results.size() == 1){
          const auto& gff_record = results.front();
          // Note: A ugly way to get gene name in attrs of gff record.
          const auto gene_name = gff_record.attrs.substr(3);
          if (expr_mat_b.contains(gene_name)) expr_mat_b[gene_name] += 1;
          else expr_mat_b[gene_name] = 1;
        }
      }
    );
  }

  // Verify expression matrix of sample A.
  CHECK(expr_mat_a["gene0"] == 20);
  CHECK(expr_mat_a["gene1"] == 10);
  CHECK(expr_mat_a["gene2"] == 20);
  CHECK(expr_mat_a["gene3"] == 10);
  CHECK(expr_mat_a["gene4"] == 20);

  CHECK(expr_mat_a["gene5"] == 0);
  CHECK(expr_mat_a["gene6"] == 0);
  CHECK(expr_mat_a["gene7"] == 0);
  CHECK(expr_mat_a["gene8"] == 0);
  CHECK(expr_mat_a["gene9"] == 0);

  // Verify expression matrix of sample B.
  CHECK(expr_mat_b["gene0"] == 20);
  CHECK(expr_mat_b["gene1"] == 10);
  CHECK(expr_mat_b["gene2"] == 20);
  CHECK(expr_mat_b["gene3"] == 10);
  CHECK(expr_mat_b["gene4"] == 20);

  CHECK(expr_mat_b["gene5"] == 0);
  CHECK(expr_mat_b["gene6"] == 0);
  CHECK(expr_mat_b["gene7"] == 0);
  CHECK(expr_mat_b["gene8"] == 0);
  CHECK(expr_mat_b["gene9"] == 0);

  // todo: sample b

  // todo: Output annotation results to BED file.
}

TEST_CASE("Tailor parallel") {
  constexpr auto READ_COUNT = 1000000;
  constexpr size_t NUM_COMPUTING_THREAD = 20;
  constexpr size_t NUM_IO_THREAD = 2;

  // Prepare index.
  auto index = bio::Index{};
  {
    auto ifs = std::ifstream{"parallel_ref.idx"};
    index.load(ifs);
  }
  auto rc_index = bio::Index{};
  {
    auto ifs = std::ifstream{"parallel_rc_ref.idx"};
    rc_index.load(ifs);
  }

  // Init tailor.
  auto tailor = bio::Tailor{index, rc_index};
  tailor.allow_seed_mismatch = true;

  // Load features into annotator.
  auto genes = bio::Annotator<bio::GffRecord>{};
  {
    auto ifs = std::ifstream{"parallel.gff"};
    auto record = bio::GffRecord{};
    while (ifs >> record) {
      // Note: Add prefix "chr" to seqname.
      if (record.seqid[0] != 'c')
        record.seqid = "chr" + record.seqid;
      genes.insert(record);
    }
  }
  genes.index(); // Note: !!!!!Important remember to index.

  // Pipeline
  {
    auto sample_list = std::vector<std::string>{"parallel_1", "parallel_2"};
    auto answer = std::map<std::string, std::vector<int>>{
      {sample_list[0], { 9,  4, 3, 4, 7, 0,  8, 5, 9, 0}},
      {sample_list[1], {16, 10, 7, 6, 4, 1, 10, 9, 7, 0}}
    };

    auto computing_task = [ // task for comp thread
      &tailor = std::as_const(tailor),
      &genes = std::as_const(genes)
    ] (std::vector<bio::FastqRecord<false>> reads, std::string sample_name, IOQueue& ioqueue) {
      // Sample alignment.
      auto aligner = [&tailor](const auto& record){
        return tailor.search(record);
      };

      // Push aligment results to IO thread
      auto sam_writer = [&ioqueue](auto&& aln){
        ioqueue.push_job(aln.first);
        return std::forward<decltype(aln)>(aln);
      };

      // Matching gene interval
      auto annotator = [&genes = std::as_const(genes)](auto&& aln){
        bool hit = (aln.first.hits.size() == 1);
        // Discard multi-mapping reads
        return hit ? genes.find(aln.first.hits.front().intv) : std::vector<bio::GffRecord>();
      };
      

      // Main pipeline
      auto batch_expr_mat = std::unordered_map<std::string, std::uint32_t>{};
      auto range_of_results = std::views::all(reads)
        | ranges::views::transform(aligner)
        | ranges::views::transform(sam_writer)
        | ranges::views::transform(annotator);
      ranges::for_each(
        range_of_results,
        [&batch_expr_mat](const auto& results){
          // Discard ambiguous annotation.
          if(results.size() == 1){
            const auto& gff_record = results.front();
            // Note: A ugly way to get gene name in attrs of gff record.
            const auto gene_name = gff_record.attrs.substr(3);
            if (batch_expr_mat.contains(gene_name)) batch_expr_mat[gene_name] += 1;
            else batch_expr_mat[gene_name] = 1;
          }
        }
      );
      return std::make_tuple(sample_name, batch_expr_mat);
    };

    auto io_task = [](std::string output_filename, IOQueue& io_queue, int sleep_time = 200){ // task for io thread
      auto ofs = std::ofstream{output_filename};
      auto local_queue = std::list<bio::Alignment>{};
      while(!io_queue.all_completed()){
        if(io_queue.queue_empty()){
          // no job has come, just sleeping
          std::this_thread::sleep_for(std::chrono::milliseconds(sleep_time));
        } else {
          // jobs stacking to the sky, should wake up to work!
          io_queue.load_job(local_queue);
        }
        // Do jobs
        while(!local_queue.empty()){
          for(auto& sam_record: aln_to_sam_list(local_queue.front()))
            ofs << sam_record << '\n';
          local_queue.pop_front();
        }
      }
      ofs.close();
      return true;
    };

    // main loop
    auto io_thread = bio::make_threadpool(NUM_IO_THREAD);
    for(auto thr_num: {1, 4, 8, 16, 24}){
      // auto batch_size = 50000;
      auto batch_size = (READ_COUNT * 2) / thr_num;
      auto total_task = 0;
      auto start = std::chrono::high_resolution_clock::now();
      
      // create IOQueue map and start IO thread
      auto io_map = std::map<std::string, IOQueue>{};
      auto io_results = std::vector<std::future<decltype(io_task("", io_map[sample_list[0]]))>>();
      for(auto& sample_name: sample_list){
        io_results.emplace_back(io_thread.submit(
          io_task,
          sample_name + ".sam",
          std::ref(io_map[sample_name]) // implicitly construct an IOQueue lock inplace
        ).second);
      }

      auto computing_thread = bio::make_threadpool(thr_num);
      auto results = std::vector<std::future<decltype(computing_task({}, "", io_map[sample_list[0]]))>>();
      for(auto& sample_name: sample_list){
        auto ifs = std::ifstream{sample_name + ".fq"};
        auto record = bio::FastqRecord<false>{};
        int count = 0;
        auto reads = std::vector<bio::FastqRecord<false>>();
        while(ifs >> record){
          reads.emplace_back(record);
          ++count;
          if(reads.size() == batch_size){
            results.emplace_back(computing_thread.submit(
              computing_task,
              reads,
              sample_name,
              std::ref(io_map[sample_name])
            ).second);
            reads.clear();
            ++total_task;
          }
        }
        if(reads.size() != 0) {
          results.emplace_back(computing_thread.submit(
            computing_task,
            reads,
            sample_name,
            std::ref(io_map[sample_name])
          ).second);
          ++total_task;
        }
      }
      // end computing task and collect results
      auto expr_mat = std::unordered_map<std::string, std::unordered_map<std::string, std::uint32_t>>{};
      auto finished_task = 0;
      while(finished_task != total_task){
        for(auto&& res_future: results){
          if(res_future.valid()){
            if(res_future.wait_for(std::chrono::microseconds(1)) == std::future_status::ready){
              auto&& [sample_name, batch_expr_mat] = res_future.get();
              if(!expr_mat.contains(sample_name))
                expr_mat[sample_name] = decltype(batch_expr_mat){};
              auto& sam_expr = expr_mat[sample_name];
              for(auto& [gene_name, expr]: batch_expr_mat){
                if (sam_expr.contains(gene_name)) sam_expr[gene_name] += expr;
                else sam_expr[gene_name] = expr;
              }
              ++finished_task;
            }
          }
        }
      }
      // end io_thread's job
      for(auto& [key, io_queue]: io_map)
        io_queue.finish_thread();
      finished_task = 0;
      while(finished_task != sample_list.size()){
        for(auto&& io_res_future: io_results){
          if(io_res_future.valid()){
            if(io_res_future.wait_for(std::chrono::microseconds(1)) == std::future_status::ready){
              bool io_finished = io_res_future.get();
              REQUIRE(io_finished == true);
              ++finished_task;
            }
          }
        }
      }
      auto end = std::chrono::high_resolution_clock::now();
      auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
      // std::cout << "Time for thr_num#" << thr_num << ": " << dur.count() << " milli s.\n";
      for(auto& sample: sample_list){
        auto& sample_expr_mat = expr_mat[sample];
        auto& ans = answer[sample];
        for(auto i = 0; i < 10; ++i){
          CHECK(sample_expr_mat["gene" + std::to_string(i)] == ans[i]);
        }          
      }
    }
  }
}