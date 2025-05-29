#include <biovoltron/algo/align/tailor/tailor.hpp>
#include <biovoltron/algo/annotate/annotator.hpp>
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
  auto expr_mat_a = std::unordered_map<std::string, std::uint32_t>{};
  auto expr_mat_b = std::unordered_map<std::string, std::uint32_t>{};
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
        aln_to_sam_list(aln),
        [&ofs](auto&& sam_record){
          ofs << sam_record << '\n';
        }
      );
      return std::forward<decltype(aln)>(aln);
    };

    // Matching gene interval
    auto annotator = [&genes](auto&& aln){
      bool hit = aln.hits.size() == 1;
      // Discard multi-mapping reads
      return hit ? genes.find(aln.hits.front().intv) : std::vector<bio::GffRecord>();
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
