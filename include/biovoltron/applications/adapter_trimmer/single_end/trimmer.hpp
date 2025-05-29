#pragma once

#include <ranges>
#include <vector>
#include <concepts>
#include <fstream>
#include <filesystem>
#include <spdlog/spdlog.h>

#include <biovoltron/file_io/all.hpp>
#include <biovoltron/utility/istring.hpp>
#include <biovoltron/algo/align/tailor/tailor.hpp>
#include <biovoltron/algo/align/tailor/index.hpp>
#include <biovoltron/algo/assemble/assembler.hpp>
#include <biovoltron/applications/adapter_trimmer/single_end/skewer/main.hpp>

using namespace std::literals;

namespace biovoltron {
namespace single {

template<class R>
class AdapterTrimmer {
private:

  using SeqView = std::conditional_t<
    R::encoded,
    istring_view,
    std::string_view
  >;
  using SeqType = std::conditional_t<R::encoded, istring, std::string>;

  class Parameter {
  public:
    /**
     * TODO: add istring constructor for std::string and std::string_view
     */
    std::string DEFAULT_ADAPTER1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"s;
    
    /**
     * the total read amount that use to detect possible adapter fragments
     */
    const std::size_t DETECT_READS = 10000;

    /**
     * The minimum length of assembled adapter
     */
    std::size_t MIN_ADAPTER_LEN = 5;

    /**
     * the minimum length for trimmed read
     */
    std::size_t MIN_READ_LEN = 0;
  } param;

  AdapterAssembler assembler;

  auto read_file(std::fstream& fin) -> std::vector<R> {
    auto read = R{};
    auto reads = std::vector<R>{};
    while (fin >> read) {
      reads.emplace_back(std::move(read));
    }
    return reads;
  }

  auto print_param(bool sensitive, std::size_t thread_num, std::size_t read_size) {
    spdlog::info("Parameter for single-end adapter trimming");
    spdlog::info("Sensitive mode: {}", sensitive);
    spdlog::info("Run with {} threads", thread_num);
    spdlog::info("Read size = {}", read_size);
    spdlog::info("Minimum read length after trimming = {}", param.MIN_READ_LEN);
    spdlog::info("Minimum assembled adapter len = {}", param.MIN_ADAPTER_LEN);
  }

  auto detect_adapter(const Tailor<>& tailor,
                      const std::vector<R>& reads,
                      const bool sensitive) {

    /* collect tails by calling Tailor */
    /* TODO use SeqView here, but we should change to SeqType */
    std::vector<SeqType> tails;
    for (int i = 0; i < (int) reads.size() && (int) tails.size() < 3000; i++) {
      auto align_res = tailor.search(reads[i]);
      if (align_res.tail_pos >= 0) {
        auto s = align_res.seq.substr(align_res.tail_pos);
        if (!s.empty()) {
          tails.emplace_back(s);
        }
      }
    }
    std::vector<SeqView> tail_views;
    std::ranges::transform(tails, std::back_inserter(tail_views), [](auto &s) {
      return SeqView(s);
    });
    auto adapter = assembler.assemble(tail_views, sensitive);
    SPDLOG_DEBUG("tail size = {}", tails.size());
    
    if (adapter.empty()) {
      spdlog::info("Unable to detect adapter, use default adapter");
      adapter = param.DEFAULT_ADAPTER1;
    }
    spdlog::info("Detected adpater = {}", adapter);
    return adapter;
  }

  auto call_skewer(const std::filesystem::path &reads_path,
                   const SeqView adapter,
                   const std::filesystem::path &output_path,
                   const std::size_t thread_num,
                   const bool sensitive) {

    spdlog::info("Call skewer for adapter trimming");

    auto skewer_out = std::filesystem::temp_directory_path() / "skewer_output";
    std::vector<std::string> skewer_argv;
    skewer_argv.emplace_back("skewer");   // main program
    skewer_argv.emplace_back(reads_path.string());
    skewer_argv.emplace_back("-o");
    skewer_argv.emplace_back(skewer_out.string());
    skewer_argv.emplace_back("-l");       /* minimum adapter length */
    skewer_argv.emplace_back(std::to_string(param.MIN_ADAPTER_LEN));
    skewer_argv.emplace_back("-t");       /* thread number */
    skewer_argv.emplace_back(std::to_string(thread_num));

    if (sensitive) {
      skewer_argv.emplace_back("-r");     /* TODO: survey what is this */
      skewer_argv.emplace_back("0.2");
    }
    if (adapter.size()) {
      skewer_argv.emplace_back("-x");
      if constexpr (R::encoded) {
        skewer_argv.emplace_back(Codec::to_string(adapter));
      } else {
        skewer_argv.emplace_back(adapter);
      }
    }

    auto skewer_argc = skewer_argv.size();
    SPDLOG_DEBUG("Called skewer: argc = {}", skewer_argc);
    {
      std::string cmd;
      std::ranges::for_each(skewer_argv, [&](auto& s) { cmd += s + ' '; });
      SPDLOG_DEBUG("skewer cmd = {}", cmd);
    }
    std::vector<const char*> real_argv(skewer_argc);
    std::ranges::for_each(std::views::iota(0u, skewer_argc), [&](auto idx) {
      real_argv[idx] = skewer_argv[idx].data();
    });
    skewer::main(skewer_argc, real_argv.data());
    skewer_out += "-trimmed.fastq";
    assert(std::filesystem::exists(skewer_out));
    std::filesystem::copy(skewer_out, output_path, 
                          std::filesystem::copy_options::overwrite_existing);
    std::filesystem::remove(skewer_out);
  }

  // In Origin EARRINGS, we have this function but never been used.
  // 
  // 
  // auto estimate_umi(const std::vector<SeqView>& tails, const SeqView adapter) {
  //   auto subseq = adapter.substr(0, 5);
  //   auto umi_sz_cnt = std::unordered_map<std::size_t, std::size_t>{};
  //   for (auto& s : tails) {
  //     for (auto i = 0u; i + subseq.size() < s.size(); i++) {
  //       if (s.substr(i, subseq.size()) == subseq) {
  //         umi_sz_cnt += 1;
  //         break;
  //       }
  //     }      
  //   }
  //   auto max_len = 0u;
  //   for (auto max_v = 0u; auto [p, cnt] : umi_sz_cnt) {
  //     if (cnt > max_v) {
  //       max_len = p;
  //       cnt = max_v;
  //     }
  //   }
  //   return max_len;
  // }

public:

  /**
   * @brief Set the default adapter, which will be used when auto-detect adapter
   * failed
   * 
   * @param adapter the adapter sequence
   * @return None
   */
  auto set_default_adapter(const std::string& adapter) {
    param.DEFAULT_ADAPTER1 = adapter;
  }

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
   * @brief Set the minimum read length after adapter trimming.
   *
   * @param len Minimum read length that the trimmed read have.
   * @return None
   */
  auto set_min_read_len(const std::size_t len) {
    param.MIN_READ_LEN = len;
  }
  

  auto trim(const Tailor<>& tailor,
            const std::vector<R>& reads,
            int thread_num = std::thread::hardware_concurrency(),
            const bool sensitive = false) {

    if (thread_num <= 0) {
      thread_num = 1;
    } else {
      thread_num = std::min(
        static_cast<int>(std::thread::hardware_concurrency()),
        thread_num
      );
    }
    const auto reads_size = reads.size();

    print_param(sensitive, thread_num, reads_size);
    auto adapter = detect_adapter(tailor, reads, sensitive);
    auto reads_path = std::filesystem::temp_directory_path() / "reads";
    {
      std::ofstream fout(reads_path);
      for (auto &r : reads) {
        fout << r << std::endl;
      }
    }
    auto output_path = 
      call_skewer(reads_path, adapter, output_path, thread_num, sensitive);

    /* Read reads in output path */
  }

  /**
   * @brief single-end adapter trimmer
   * @param tailor 
   * @param reads_path the read file path
   * @param output_path the output file path
   * @param thread_num the thread amount, default is
   * `std::thread::hardware_concurrenty()`
   * @param sensitive trim read in sensitive mode
   */
  void trim(const Tailor<>& tailor,
            const std::filesystem::path reads_path,
            const std::filesystem::path output_path,
            int thread_num = std::thread::hardware_concurrency(),
            const bool sensitive = false) {
  
    auto open_file = [](const std::filesystem::path& path, const bool is_input) {
      /* check whether is exists */
      if (!std::filesystem::exists(path)) {
        SPDLOG_ERROR("{} doesn't exists", path.string());
        std::exit(1);
      }
      std::fstream f(path, is_input ? std::ios::in : std::ios::out);
      if (!f.is_open()) {
        SPDLOG_ERROR("Cannot open {}", path.string());
        std::exit(-1);
      }
      return f;
    };

    auto read_fin = open_file(reads_path, true);
    auto reads = read_file(read_fin);
    print_param(sensitive, thread_num, reads.size());
    auto adapter = detect_adapter(tailor, reads, sensitive);
    call_skewer(reads_path, adapter, output_path, thread_num, sensitive);
  }

};

} // namespce single
} // namespace biovoltron