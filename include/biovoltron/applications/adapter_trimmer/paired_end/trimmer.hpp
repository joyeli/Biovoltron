#pragma once

#include <thread>
#include <ranges>
#include <random>
#include <filesystem>
#include <sstream>
#include <fstream>
#include <map>
#include <set>

#include <boost/asio.hpp>
#include <spdlog/spdlog.h>
#include <omp.h>

#include <biovoltron/applications/adapter_trimmer/detail/simd.hpp>
#include <biovoltron/algo/assemble/assembler.hpp>
#include <biovoltron/file_io/all.hpp>
#include <biovoltron/utility/istring.hpp>
#include <biovoltron/container/xbit_vector.hpp>


using namespace std::literals;

namespace biovoltron {

/**
 * @ingroup applications
 * @brief The paired-end reads adapter trimmer, here is the implementation from 
 * EARRINGS, which optimized by multithreading and SIMD intrinsic functions.
 * @todo single-end reads adapter trimmer -> wait tailor
 * ```cpp
 * #include <biovoltron/applications/adapter_trimmer/trimmer.hpp>
 * #include <filesystem>
 * #include <fstream>
 * 
 * int main() {
 *  TODO: usage
 *   std:: 
 * }
 * ```
 * 
 * @tparam R the type of reads, we only take FastaRecord and FastqRecord now,
 * will support BamRecord and SamRecord in future.
 * @ref https://academic.oup.com/bioinformatics/article/37/13/1846/6103563?login=true
 */
namespace paired {

template<class R>
class AdapterTrimmer {

  // TODO: modify if add SamRecord support
  using SeqType = std::conditional_t<R::encoded, istring, std::string>;
  using SeqView = std::conditional_t<
    R::encoded,
    istring_view,
    std::string_view
  >;

  class Parameter {
  public:
    /**
     * TODO: add istring constructor for std::string and std::string_view
     **/
    std::string DEFAULT_ADAPTER1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"s;
    std::string DEFAULT_ADAPTER2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA"s;

    /**
     * [RC check]
     * The compare length use to check the similarity between reverse 
     * complement of forward reads and reverse reads, and vice versa
     */
    std::size_t RC_CHECK_LEN = 16;
    
    /**
     * The read size use to detect possible adapter fragments, which will use
     * to assemble adapters.
     */
    const std::size_t DETECT_READS = 10000;
    
    /**
     * [RC check]
     * the similarity threshold that considered hit if the similarity of
     * `RC_CHECK_LEN` base from 5' of R1 and R2 is greater or equal to this
     * value
     */
    double TAIL_MATCH_RATIO = 0.7;

    /**
     * [SS check]
     * the similarity threshold that considered hit if the similarity of
     * DNA insert part of R1 and R2 is greater than or equal to this value
     */
    double DNA_MATCH_RATIO = 0.9;

    /**
     * [AS check]
     * the similarity threshold that considered hit if the similarity of
     * adapter part and assembled adapter is greater than or equal to this
     * value
     */
    double ADAPTER_MATCH_RATIO = 0.8;
  
    /**
     * the minimum length of assembled adapter
     */
    std::size_t MIN_ADAPTER_LEN = 5;

    /**
     * the checked length when adapter is located on 3'
     */
    std::size_t HEAD_ADAPTER_CHECK_LEN = 16;

    /**
     * the minimum length for trimmed read
     */
    std::size_t MIN_READ_LEN = 0;

    /**
     * the size of a single asio read batch
     */
    std::size_t ASIO_BUFFER_SIZE = 8192;
  } param;

  AdapterAssembler assembler;

  auto print_param(const bool sensitive, const std::size_t thread_num, bool asio_mode = false) {
    spdlog::info("Paired-end adapter trimming");
    spdlog::info("Asio mode = {}", asio_mode);
    spdlog::info("Sensitive mode: {}", sensitive);
    spdlog::info("Run with {} threads", thread_num);
    spdlog::info("Tail match ratio = {:.3f}", param.TAIL_MATCH_RATIO);
    spdlog::info("DNA match ratio = {:.3f}", param.DNA_MATCH_RATIO);
    spdlog::info("Adapter match ratio = {:.3f}", param.ADAPTER_MATCH_RATIO);
    spdlog::info("Minimum read length after trimming = {}", param.MIN_READ_LEN);
    spdlog::info("Asio buffer size = {}", param.ASIO_BUFFER_SIZE);
    if (asio_mode) {
    #ifdef SIMDPP_ARCH_X86_AVX2 
      spdlog::info("SIMD instruction set used : AVX2");
    #else
      spdlog::info("SIMD instruction set used : SSE3");
    #endif
    }
    SPDLOG_DEBUG("Vector size= {}", detail::VECTOR_SIZE);
    SPDLOG_DEBUG("BASE in VECTOR = {}", detail::BASE_IN_VECTOR);
  }

  /**
   * @brief Given sequence and possible adapter, get all possible positions of
   * adapter appearence in the sequence. 
   * @details Calculate the similarity between all possible subsequences and 
   * adapter, which considered hit if it is over `TAIL_MATCH_RATIO`.
   * 
   * @param seq the sequence
   * @param adapter possible adapter
   * @return std::vector<std::size_t> 
   */
  auto find_possible_adapter_pos(SeqView seq,
                                 SeqView adapter) const noexcept 
                                 -> std::vector<std::size_t> {
    auto possible_pos = std::vector<std::size_t> {};

    const auto read_size = std::min(detail::BASE_IN_VECTOR, seq.size());
    const auto cmp_len = std::min(param.RC_CHECK_LEN, read_size);
    const auto subseq = seq.substr(0, read_size);
    auto seq_vec = detail::make_simd_vector(subseq, 1)[0];
    const auto adapter_vec = detail::make_simd_vector(adapter, 1)[0];

    /* interate all possible subsequence */
    for (auto i = 0U; i <= read_size - cmp_len; i++) {
      auto simd_sim = detail::cal_similarity(seq_vec, adapter_vec, cmp_len);
      if (simd_sim >= param.TAIL_MATCH_RATIO) {
        possible_pos.emplace_back(i + cmp_len);
      }
      seq_vec = detail::shift_left(seq_vec);
    }
    /**
     * The base number of seq is greater than detail::BASE_IN_VECTOR, which means we
     * will have more than one SIMD_VECTOR to save to sequence. In order to
     * shift more than one SIMD_VECTOR, we use recursive function to get
     * possible position after detail::BASE_IN_VECTOR - RC_CHECK_LEN
     */
    if (seq.size() > detail::BASE_IN_VECTOR) {
      auto base_pos = detail::BASE_IN_VECTOR - param.RC_CHECK_LEN + 1;
      auto&& r = find_possible_adapter_pos(seq.substr(base_pos), adapter);
      std::ranges::transform(r, std::back_inserter(possible_pos), 
                             [&](const auto& idx) {
                               return idx + base_pos;
                             });
    }
    return possible_pos;
  }

  /**
   * @brief After we have intersect positions from both reads, we want to find
   * the closest position to 5' that seperate the sequence into DNA part and
   * adapter part.
   * @details We perform Sequence similarity(SS) check, which calculated the 
   * similarity between DNA part from both reads after determine the seperation
   * position, we considered a hit if it is greater or equal to
   * `DNA_MATCH_RATIO`.
   * 
   * @param forward_seq forward sequence
   * @param reverse_seq reverse sequence
   * @param intersect_pos intersect of possible adapter positions
   * @return `std::size_t`
   */
  auto get_adapter_pos(const SeqView forward_seq, const SeqView reverse_seq,
                       const std::vector<std::size_t>& intersect_pos) 
                       const noexcept
                       -> std::size_t {

    auto forward_rc = Codec::rev_comp(forward_seq);

    /* we only need to make simd vector of reverse read for once */
    auto rev_vec = detail::make_simd_vector(reverse_seq);

    for (auto fwd_rc_view = SeqView(forward_rc);
         auto pos : intersect_pos | std::views::reverse) {

      auto fwd_rc_subseq = fwd_rc_view.substr(fwd_rc_view.size() - pos);
      auto fwd_rc_subseq_vec = detail::make_simd_vector(fwd_rc_subseq);
      auto match = 0;
      for (auto i = 0u; i < fwd_rc_subseq_vec.size(); i++) {
        auto cmp_len = std::min(
          detail::BASE_IN_VECTOR,
          pos - i * detail::BASE_IN_VECTOR
        );
        match += detail::cal_match(fwd_rc_subseq_vec[i], rev_vec[i], cmp_len);
      }
      if (static_cast<double>(match) / pos >= param.DNA_MATCH_RATIO) {
        return pos;
      }
    }
    return static_cast<std::size_t>(-1);
  }

  /**
   * @brief find proper trimming position for each read
   * @details After we get assembled adatper and intersecting positions between 
   * both reads. For every intersection, we will first check whether the adapter 
   * part pass AS check,then check whether the DNA part pass SS check. If so, 
   * the position which closest to 5' is considered as trimming position.
   * 
   * If none of intersecting position meet the condition above, we will check 
   * whether the adapter is located on 3'. By checking the 
   * 
   * 
   *                   pos
   * 3' ----------------|---------------------- 5'
   *        DNA part          adapter part
   * 
   * @param forward_seq forward sequence
   * @param reverse_seq reverse sequence
   * @param intersect_pos intersecting position of possible adapter
   * @param forward_adapter forward assembled adapter
   * @param reverse_adapter reverse assembled adapter
   * @return `std::size_t`
   */
  auto get_trim_pos(const SeqView forward_seq, const SeqView reverse_seq,
                    const std::vector<std::size_t>& intersect_pos,
                    const SeqView forward_adapter,
                    const SeqView reverse_adapter) 
                    const noexcept -> std::size_t {

    auto forward_rc = Codec::rev_comp(forward_seq);
    
    /* we only need to make simd vector of reverse read for once */
    auto rev_vec = detail::make_simd_vector(reverse_seq);
    auto fwd_adapter_vec = detail::make_simd_vector(forward_adapter, 1)[0];
    auto rev_adapter_vec = detail::make_simd_vector(reverse_adapter, 1)[0];

    /* like `get_adapter_pos()`, we perform SS check here */
    for (auto fwd_rc_view = SeqView(forward_rc);
         auto pos : intersect_pos | std::views::reverse) {
      if (pos == forward_seq.size()) {
        continue;
      }
      auto fwd_rc_subseq = fwd_rc_view.substr(forward_rc.size() - pos);
      auto fwd_rc_subseq_vec = detail::make_simd_vector(fwd_rc_subseq);

      auto match = 0;
      for (auto i = 0u; i < fwd_rc_subseq_vec.size(); i++) {
        auto cmp_len = std::min(
          detail::BASE_IN_VECTOR,
          pos - i * detail::BASE_IN_VECTOR
        );
        match += detail::cal_match(fwd_rc_subseq_vec[i], rev_vec[i], cmp_len);
      }

      /* pass the SS check, perform AS check */
      if (static_cast<double>(match) / pos >= param.DNA_MATCH_RATIO) {
        auto fwd_subseq = forward_seq.substr(pos);
        auto rev_subseq = reverse_seq.substr(pos);

        auto fwd_subseq_vec = detail::make_simd_vector(fwd_subseq, 1)[0];
        auto rev_subseq_vec = detail::make_simd_vector(rev_subseq, 1)[0];
        auto len = std::min(fwd_subseq.size(), forward_adapter.size());
        
        /**
         * if the similarity between adapter part and assembled adapter is
         * greater or equal to `ADAPTER_MATCH_RATIO`, then this position is a 
         * proper position for trimming.
         */
        if (detail::cal_similarity(fwd_subseq_vec, fwd_adapter_vec, len)
            >= param.ADAPTER_MATCH_RATIO ||
            detail::cal_similarity(rev_subseq_vec, rev_adapter_vec, len)
            >= param.ADAPTER_MATCH_RATIO) {
          return pos;
        }
      }
    }

    /**
     * No proper trimming position is detected in reads, check adapter is
     * located at the head of reads
     */
    if (intersect_pos.empty()) {
      auto find_head_adapter_pos = [&](const SeqView seq,
                                       const detail::SIMD_VECTOR& adapter_vec,
                                       const std::size_t adapter_size) {
        auto seq_vec = detail::make_simd_vector(seq, 1)[0];
        auto cmp_len = std::min(param.HEAD_ADAPTER_CHECK_LEN, adapter_size);
        for (auto pos : std::views::iota(0UL, param.HEAD_ADAPTER_CHECK_LEN)) {
          cmp_len = std::min(seq.size() - pos, cmp_len);
          if (cmp_len < param.MIN_ADAPTER_LEN) {
            break;
          }
          auto similarity = detail::cal_similarity(
            seq_vec, adapter_vec, cmp_len
          );
          if (similarity >= param.ADAPTER_MATCH_RATIO) {
            return pos;
          }
          seq_vec = detail::shift_left(seq_vec);
        }
        return static_cast<std::size_t>(-1);
      };

      // ? figure out the origin implementation look like this
      // from 32 to min_adapter_len, check the similarity of both sequence, 
      // check that whether adapter is at the head, then trim it.
      {
        auto forward_subseq = forward_seq.substr(
            0, 2 * param.HEAD_ADAPTER_CHECK_LEN);
        auto pos = find_head_adapter_pos(forward_subseq, fwd_adapter_vec,
                                         forward_adapter.size());
        if (pos != static_cast<std::size_t>(-1)) {
          return pos;
        }
      }
      {
        auto reverse_subseq = reverse_seq.substr(
            0, 2 * param.HEAD_ADAPTER_CHECK_LEN);
        auto pos = find_head_adapter_pos(reverse_subseq, rev_adapter_vec,
                                         reverse_adapter.size());
        if (pos != static_cast<std::size_t>(-1)) {
          return pos;
        }
      }
    }

    return forward_seq.size();
  }

  /**
   * @brief preprocess for paired-end reads, resize both reads to same length,
   * and transform `N` to randomly "ATCG"
   * 
   * @param fwd_read forward read
   * @param rev_read reverse read
   * @return None
   */
  auto preprocess(R& fwd_read, R& rev_read) const noexcept -> void {
    auto sz = std::min(fwd_read.seq.size(), rev_read.seq.size());
    fwd_read.seq.resize(sz);
    rev_read.seq.resize(sz);
    if constexpr (std::same_as<R, FastqRecord<R::encoded>>) {
      fwd_read.qual.resize(sz);
      fwd_read.qual.resize(sz);
    }
    static std::mt19937 rng(
      std::chrono::steady_clock::now().time_since_epoch().count()
    );
    auto check_for_N = [](auto c) {
      if constexpr (R::encoded) {
        if (c > 3) {
          c = rng() & 3;
        }
      } else {
        if (c == 'N') {
          c = Codec::to_char(rng() & 3);
        }
      }
      return c;
    };
    std::ranges::transform(fwd_read.seq, fwd_read.seq.begin(), check_for_N);
    std::ranges::transform(rev_read.seq, rev_read.seq.begin(), check_for_N);
  }


  /**
   * @brief get the possible adapter positions by calcuate the similarity 
   * between reverse complement of reverse read and forward read, and vice versa.
   * We will return the intersect between two set of positions.
   * @param fwd_seq forward sequence
   * @param rev_seq reverse sequence
   * @return std::vector<std::size_t>
   */
  auto get_possible_intersect(SeqView fwd_seq,
                              SeqView rev_seq) const noexcept 
                              -> std::vector<std::size_t> {
    auto seq_pair = std::array { fwd_seq, rev_seq };
    auto possible_pos = std::array<std::vector<std::size_t>, 2>{};

    for (auto s = 0; s <= 1; s++) {
      auto rc = Codec::rev_comp(seq_pair[s].substr(0, param.RC_CHECK_LEN));
      possible_pos[s] = find_possible_adapter_pos(seq_pair[s ^ 1], rc);
    }

    auto intersect_pos = std::vector<std::size_t> {};
    std::ranges::set_intersection(possible_pos[0],
                                  possible_pos[1],
                                  std::back_inserter(intersect_pos));
    return intersect_pos;
  };
  
  /**
   * @brief detect adatper from the first `DETECT_READS` reads of given files.
   * 
   * @param forward_reads 
   * @param reverse_reads 
   * @param sensitive 
   * @return std::pair<SeqType, SeqType>
   */
  auto detect_adapter(const std::vector<R>& forward_reads,
                      const std::vector<R>& reverse_reads,
                      bool sensitive) noexcept 
                      -> std::pair<SeqType, SeqType> {

    auto forward_tails = std::vector<SeqView> {};
    auto reverse_tails = std::vector<SeqView> {};
    auto sz = forward_reads.size();

    spdlog::info("Collect tails from forward reads and reverse reads...");
    for (auto i = 0u; i < sz; i++) {
      auto fwd_seq = SeqView(forward_reads[i].seq);
      auto rev_seq = SeqView(reverse_reads[i].seq);
      auto intersect_pos = get_possible_intersect(fwd_seq, rev_seq);

      /* No intersect positions, skip this reads */
      if (intersect_pos.size() == 0) {
        continue;
      }

      auto pos = get_adapter_pos(fwd_seq, rev_seq, intersect_pos);

      /* cannot find proper adapter position */
      if (pos == static_cast<std::size_t>(-1) || pos == fwd_seq.size()) {
        continue;
      }

      forward_tails.push_back(fwd_seq.substr(pos));
      reverse_tails.push_back(rev_seq.substr(pos));
    }
    spdlog::info("Done. Start detect adapter from collected tails");
    SPDLOG_DEBUG("forward tails size = {}", forward_tails.size());
    SPDLOG_DEBUG("reverse tails size = {}", reverse_tails.size());

    spdlog::info("Detect forward adapter...");
    auto forward_adapter = assembler.assemble(forward_tails, sensitive);
    spdlog::info("Detect reverse adapter...");
    auto reverse_adapter = assembler.assemble(reverse_tails, sensitive);

    /* unable to assemble adapter, use default adapter */
    if (forward_adapter.empty() || reverse_adapter.empty()) {
      spdlog::info("Unable to detect adapter in reads");
      spdlog::info("Use default adapter");
      if constexpr (R::encoded) {
        forward_adapter = Codec::to_istring(param.DEFAULT_ADAPTER1);
        reverse_adapter = Codec::to_istring(param.DEFAULT_ADAPTER2);
      } else {
        forward_adapter = param.DEFAULT_ADAPTER1;
        reverse_adapter = param.DEFAULT_ADAPTER2;
      }
    }

    /* trim the adapter to same length */
    if (forward_adapter.size() > reverse_adapter.size()) {
      forward_adapter.resize(reverse_adapter.size());
    } else if (forward_adapter.size() < reverse_adapter.size()) {
      reverse_adapter.resize(forward_adapter.size());
    }

    spdlog::info("The forward adapter = {}", forward_adapter);
    spdlog::info("The reverse adapter = {}", reverse_adapter);

    return std::make_pair(forward_adapter, reverse_adapter);
  }

  /**
   * @brief After we get assembled adatper, trim reads on the postion by calling
   * `get_trim_pos`
   * 
   * @param forward_read forward read
   * @param reverse_read reverse read
   * @param forward_adapter forward assembled adapter
   * @param reverse_adapter reverse assembled adapter
   * @return void
   */
  auto trim_read(R& forward_read, R& reverse_read,
                 const SeqType& forward_adapter, const SeqType& reverse_adapter) 
                 const noexcept
                 -> void {
    auto fwd_seq_view = SeqView(forward_read.seq);
    auto rev_seq_view = SeqView(reverse_read.seq);

    auto intersect_pos = get_possible_intersect(fwd_seq_view, rev_seq_view);

    auto trim_pos = get_trim_pos(fwd_seq_view, rev_seq_view,
                                 intersect_pos,
                                 forward_adapter, reverse_adapter);

    if (trim_pos < param.MIN_READ_LEN) {
      trim_pos = 0;
    }
    forward_read.seq.resize(trim_pos);
    reverse_read.seq.resize(trim_pos);
    if constexpr (std::same_as<R, FastqRecord<R::encoded>>) {
      forward_read.qual.resize(trim_pos);
      reverse_read.qual.resize(trim_pos);
    }
  }

  /**
   * @brief trim reads in asynchronous mode
   * @details
   * 
   * @param forward_fin std::ifstream of file contains forward reads
   * @param reverse_fin std::ifstream of file contains reverse reads
   * @param forward_fout std::ofstream of file store trimmed forward reads
   * @param reverse_fout std::ofstream of file store trimmed reverse reads
   * @param thread_num the number of threads
   * @param sensitive whether perform trimming in sensitive mode
   * @return `void`
   */
  auto asio_trim(std::ifstream& forward_fin, std::ifstream& reverse_fin,
                 std::ofstream& forward_fout, std::ofstream& reverse_fout,
                 int thread_num, const bool sensitive) {

    using BufType = std::vector<R>;

    auto tmp_dir = std::filesystem::temp_directory_path();
    auto fwd_tmp_dir = tmp_dir / "fwd";
    auto rev_tmp_dir = tmp_dir / "rev";
    std::filesystem::remove_all(fwd_tmp_dir);
    std::filesystem::remove_all(rev_tmp_dir);
    if (!std::filesystem::create_directory(fwd_tmp_dir) ||
        !std::filesystem::create_directory(rev_tmp_dir)) {
      SPDLOG_ERROR("Cannot create temporary directory {} and {}", fwd_tmp_dir.string(), rev_tmp_dir.string());
      std::exit(1);
    }
    auto ios = boost::asio::io_context {};
    auto guard = boost::asio::make_work_guard(ios);
    auto total_batch = -1LL;
    auto writed_batch = std::set<std::size_t> {};
    auto writed_batch_mtx = std::mutex {};

    std::vector<std::jthread> threads(thread_num - 1);
    for (auto i = 0; i < thread_num - 1; i++) {
      threads.emplace_back(std::jthread([&ios]() {
        ios.run();
      }));
    }

    SPDLOG_DEBUG("Start Trimming by using Asychronous IO");

    auto write_all_files = [&]() {
      auto batch = 0LL;
      auto final_write = [](std::filesystem::path p, std::ofstream& fout) {
        auto fin = std::ifstream(p, std::ios::binary);
        assert(fin.is_open());
        fout << fin.rdbuf() << std::flush;
        fin.close();
        std::filesystem::remove(p);
      };
      while (1) {
        while (writed_batch.contains(batch)) {
          auto fwd_path = fwd_tmp_dir / std::to_string(batch);
          auto rev_path = rev_tmp_dir / std::to_string(batch);
          // SPDLOG_DEBUG("[write all] write {}, {}", fwd_path, rev_path);
          final_write(fwd_path, forward_fout);
          final_write(rev_path, reverse_fout);
          batch++;
        }
        if (static_cast<int>(writed_batch.size()) == total_batch) {
          assert(total_batch == batch);
          // SPDLOG_DEBUG("[write all] Write files finished, batch = {}", batch);
          break;
        }
        std::this_thread::sleep_for(1ms);
      }
    };

    auto asio_read = [&](const std::size_t read_count) {
      auto fwd_buf = BufType();
      auto rev_buf = BufType();
      fwd_buf.reserve(read_count);
      rev_buf.reserve(read_count);
      bool fwd_eof = false, rev_eof = false;
      for (auto i = 0u; i < read_count; i++) {
        if (R forward_read; !(forward_fin >> forward_read)) {
          fwd_eof = true;
        } else {
          fwd_buf.emplace_back(std::move(forward_read));
        }
        if (R reverse_read; !(reverse_fin >> reverse_read)) {
          rev_eof = true;
        } else {
          rev_buf.emplace_back(std::move(reverse_read));
        }
        if (fwd_eof != rev_eof) {
          SPDLOG_ERROR("The forward reads file and reverse reads files doesn't contain same amount reads");
          exit(EXIT_FAILURE);
        } else if (fwd_eof) {
          return std::make_tuple(true, std::move(fwd_buf), std::move(rev_buf));
        }
      }
      return std::make_tuple(false, std::move(fwd_buf), std::move(rev_buf));
    };

    auto asio_write = [&](const std::size_t batch,
                          const BufType& fwd_buf, const BufType& rev_buf) {
      auto write_to_tmp = [](std::filesystem::path p, const BufType& buf) {
        auto fout = std::ofstream(p, std::ios::binary);
        if (!fout.is_open()) {
          SPDLOG_ERROR("Cannot open tmp file {}", p.string());
          exit(EXIT_FAILURE);
        }
        for (const auto& r : buf) {
          fout << r << '\n';
        }
      };

      auto fwd_path = fwd_tmp_dir / std::to_string(batch);
      auto rev_path = rev_tmp_dir / std::to_string(batch);
      write_to_tmp(fwd_path, fwd_buf);
      write_to_tmp(rev_path, rev_buf);
      auto lock = std::scoped_lock<std::mutex>(writed_batch_mtx);
      writed_batch.insert(batch);
    };

    auto&& [_, fwd_reads, rev_reads] = asio_read(param.DETECT_READS);
    auto [forward_adapter, reverse_adapter]
      = detect_adapter(fwd_reads, rev_reads, sensitive);

    auto trim_reads = [&](const std::size_t batch,
                          BufType& fwd_buf, BufType& rev_buf) {
      auto sz = fwd_buf.size();
      for (auto i = 0u; i < sz; i++) {
        trim_read(fwd_buf[i], rev_buf[i], forward_adapter, reverse_adapter);
      }
      asio_write(batch, fwd_buf, rev_buf);
    };

    auto pipeline = [&]() {
      auto batch = 0;
      auto bg = std::chrono::high_resolution_clock::now();
      while (true) {
        auto&& [eof, fwd_buf, rev_buf] = asio_read(param.ASIO_BUFFER_SIZE);
        if (fwd_buf.empty()) {
          break;
        }
        // SPDLOG_DEBUG("[Pipeline] Read batch {} finished", batch);
        boost::asio::post(ios, std::bind(trim_reads, batch,
                                         std::move(fwd_buf),
                                         std::move(rev_buf)));
        batch++;
        if (eof) {
          forward_fin.close();
          reverse_fin.close();
          break;
        }
      }
      auto ed = std::chrono::high_resolution_clock::now();
      SPDLOG_DEBUG("[Asio trim] read files takes {} ms", (ed - bg) / 1ms);
      total_batch = batch;
      
      boost::asio::post(ios, [&guard]() { guard.reset(); });
    };

    forward_fin.seekg(0);
    reverse_fin.seekg(0);

    spdlog::info("Start trimming...");
    if (thread_num <= 4) {
      pipeline();
      ios.run();
      write_all_files();
    } else {
      boost::asio::post(ios, pipeline);
      boost::asio::post(ios, write_all_files);
      ios.run();
      // pipeline();
      // ios.run();
      // write_all_files();
    }
    spdlog::info("Done. Start clear temporary directory");
    std::filesystem::remove_all(fwd_tmp_dir);
    std::filesystem::remove_all(rev_tmp_dir);
  }

public:

  /**
   * @brief Set the prune factor which used in assembing adapter, if edge
   * frequency is less than the factor in deBurjin graph, then it will not be
   * considered when expand the path.
   *
   * @param new_factor
   * @return auto
   */
  auto set_prune_factor(const double new_factor) {
    assembler.set_prune_factor(new_factor);
  }

  /**
   * @brief Set the RC check threshold. The RC check threshold is use
   * 
   * @param ratio RC check threshold, the value should between [0.0 - 1.0]
   * @return None
   */
  auto set_rc_check_threshold(const double ratio) {
    if (ratio < 0 || 1 < ratio) {
      SPDLOG_ERROR("The value of RC check ratio must inside [0 - 1], use default value {}", param.TAIL_MATCH_RATIO);
    } else {
      param.TAIL_MATCH_RATIO = ratio;
    }
  }

  /**
   * @brief Set the SS check threshold. 
   * 
   * @param ratio SS check threshold, the value should between [0.0 - 1.0]
   * @return None 
   */
  auto set_ss_check_threshold(const double ratio) noexcept {
    if (ratio < 0 || 1 < ratio) {
      SPDLOG_ERROR("The value of RC check ratio must inside [0 - 1], use default value {}", param.DNA_MATCH_RATIO);
    } else {
      param.DNA_MATCH_RATIO = ratio;
    }
  }

  /**
   * @brief Set the AS check threshold.
   * 
   * @param ratio AS check threshold, the value should between [0.0 - 1.0]
   * @return None 
   */
  auto set_as_check_threshold(const double ratio) noexcept {
    if (ratio < 0 || 1 < ratio) {
      SPDLOG_ERROR("The value of AS check ratio must inside [0 - 1], use default value {}", param.ADAPTER_MATCH_RATIO);
    } else {
      param.ADAPTER_MATCH_RATIO = ratio;
    }
  }

  /**
   * @brief Set the minimum read length after adapter trimming.
   * 
   * @param len Minimum read length that the trimmed read have.
   * @return None
   */
  auto set_min_read_len(const std::size_t len) noexcept {
    param.MIN_READ_LEN = len;
  }


  /**
   * @brief Set the default adapter1, which will be used when cannot auto-detect
   * adapter1
   * 
   * @param adapter 
   * @return None 
   */
  auto set_default_adapter1(const std::string& adapter) {
    param.DEFAULT_ADAPTER1 = adapter;
  }

  /**
   * @brief Set the default adapter2, which will be used when cannot auto-detect
   * adapter2
   * 
   * @param adapter 
   * @return None 
   */
  auto set_default_adapter2(const std::string& adapter) {
    param.DEFAULT_ADAPTER2 = adapter;
  }

  /**
   * @brief Set the asio buffer size. This value is being used in ASIO mode
   * trimming.
   * 
   * @param size The asio buffer size.
   * @return `void`
   */
  auto set_asio_buffer_size(const std::size_t size) noexcept {
    param.ASIO_BUFFER_SIZE = size;
  }

  /**
   * @brief Paired-end read adapter trimming in asychronous I/O(ASIO) mode. 
   * Compare to non-asychronous I/O mode, ASIO mode consume less memory and use 
   * less cpu usage when trimming, but it will increase the overload on 
   * system I/O.
   * 
   * @param forward_reads_path the file path contains forward reads
   * @param reverse_reads_path the file path contains reverse reads
   * @param forward_output_path the file path to store the trimmed forward reads
   * @param reverse_output_path the file path to store the trimmed reverse reads
   * @param thread_num the maximum thread used to trimming the reads
   * @param sensitive whether trim the reads in sensitive mode. The sensitive
   * mode should be used when guaranted the forward reads and reverse reads 
   * contains adapters.
   * @return None
   */
  auto trim(const std::filesystem::path forward_reads_path,
            const std::filesystem::path reverse_reads_path,
            const std::filesystem::path forward_output_path,
            const std::filesystem::path reverse_output_path,
            int thread_num = std::thread::hardware_concurrency(),
            const bool sensitive = false) {

    auto check_path_exists = [](const std::filesystem::path& path) {
      if (!std::filesystem::exists(path)) {
        SPDLOG_ERROR("{} doesn't exists", path.string());
        std::exit(1);
      }
    };

    auto open_input_file = [](const std::filesystem::path& path) {
      std::ifstream fin(std::filesystem::absolute(path));
      if (!fin.is_open()) {
        SPDLOG_ERROR("Cannot open {}", path.string());
        std::exit(1);
      }
      return fin;
    };

    auto open_output_file = [](const std::filesystem::path& path) {
      std::ofstream fout(std::filesystem::absolute(path), std::ios::binary);
      if (!fout.is_open()) {
        SPDLOG_ERROR("Cannot open {}", path.string());
        std::exit(1);
      }
      return fout;
    };

    check_path_exists(forward_reads_path);
    check_path_exists(reverse_reads_path);
    auto fwd_fin = open_input_file(forward_reads_path);
    auto rev_fin = open_input_file(reverse_reads_path);
    auto fwd_fout = open_output_file(forward_output_path);
    auto rev_fout = open_output_file(reverse_output_path);

    if (thread_num <= 0) {
      thread_num = 1;
    } else {
      thread_num = std::min(int(std::thread::hardware_concurrency()),
                            thread_num);
    }

    print_param(sensitive, thread_num, true);
    asio_trim(fwd_fin, rev_fin, fwd_fout, rev_fout, thread_num, sensitive);
  }

  /**
   * @brief paired-end adapter trimming.
   *
   * @param forward_reads `std::vector<R>` contains all forward reads
   * @param reverse_reads `std::vector<R>` contains all reverse reads
   * @param thread_num the maximum thread used to trimming the reads
   * @param sensitive whether perfrom adapter trimming in sensitive mode. The 
   * sensitive mode can be used when you guarantee there are adapter inside the
   * reads
   * @return `std::pair<std::vector<R>, std::vector<R>>` contains trimmed
   * forward reads and reverse reads
   */
  auto trim(std::vector<R> forward_reads, std::vector<R> reverse_reads,
            int thread_num = std::thread::hardware_concurrency(),
            const bool sensitive = false) {
    
    if (forward_reads.size() != reverse_reads.size()) {
      SPDLOG_ERROR("The size between forward read and reverse read is differnt");
      exit(EXIT_FAILURE);
    }
    if (thread_num <= 0) {
      thread_num = 1;
    } else {
      thread_num = std::min(
        static_cast<int>(std::thread::hardware_concurrency()),
        thread_num
      );
    }

    /* the reads count */
    const auto reads_size = forward_reads.size();
    print_param(sensitive, thread_num, false);
    spdlog::info("The read size = {}", reads_size);
    #pragma omp parallel for
    for (auto i = 0u; i < reads_size; i++) {
      preprocess(forward_reads[i], reverse_reads[i]);
    }

    auto detect_reads = std::min(reads_size, param.DETECT_READS);
    auto fwd_sub_reads = std::vector<R>(forward_reads.begin(),
                                        forward_reads.begin() + detect_reads);
    auto rev_sub_reads = std::vector<R>(reverse_reads.begin(),
                                        reverse_reads.begin() + detect_reads);
    auto [fwd_adapter, rev_adapter] = detect_adapter(fwd_sub_reads,
                                                     rev_sub_reads,
                                                     sensitive);
    
    #pragma omp parallel for
    for (auto i = 0ul; i < reads_size; i++) {
      trim_read(forward_reads[i], reverse_reads[i], fwd_adapter, rev_adapter);
    }
    return std::make_pair(forward_reads, reverse_reads);
  }
};

} // namespace paired


} // namespace biovoltron