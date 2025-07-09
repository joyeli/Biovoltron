#pragma once

#include <cstdio>

#include <spdlog/spdlog.h>

#include <intel_native/pairhmm/avx-pairhmm.h>
#include <intel_native/pairhmm/shacc_pairhmm.h>
#include <omp.h>

#include <algorithm>
#include <biovoltron/algo/align/inexact_match/pairhmm.hpp>
#include <biovoltron/file_io/sam.hpp>
#include <biovoltron/utility/haplotype/haplotype.hpp>
#include <cmath>

namespace biovoltron {

/**
 * @ingroup align
 *
 * @brief TBA
 */
struct AvxPairHMM {

 private:
  bool g_use_double;
  int g_max_threads;

  Context<float> g_ctxf;
  Context<double> g_ctxd;

  float (*g_compute_full_prob_float)(testcase* tc);
  double (*g_compute_full_prob_double)(testcase* tc);

  std::vector<std::vector<testcase>> m_testcases;

  inline auto getData(const std::vector<SamRecord<>>& readDataArray,
          const std::vector<Haplotype>& haplotypeDataArray) {
    const auto numReads = readDataArray.size();
    const auto numHaplotypes = haplotypeDataArray.size();

    auto haplotypes = std::vector<const char*>{};
    auto haplotypeLengths = std::vector<int>{};

    auto total_hap_length = 0;
    auto total_read_length = 0;

    // get haplotypes
    for (auto i = 0; i < numHaplotypes; i++) {
      const auto length = haplotypeDataArray[i].seq.length();
      haplotypes.push_back(haplotypeDataArray[i].seq.data());
      haplotypeLengths.push_back(length);
      total_hap_length += length;
    }

    // get reads and create testcases
    for (auto r = 0; r < numReads; r++) {
      const auto length = readDataArray[r].seq.length();
      const auto* reads = readDataArray[r].seq.data();
      const auto readLength = static_cast<int>(length);
      const char* insGops = readDataArray[r].insertion_gop().data();
      const char* delGops = readDataArray[r].deletion_gop().data();
      const char* gapConts = readDataArray[r].overall_gcp().data();
      const char* readQuals = readDataArray[r].qual.data();
      total_read_length += length;

      auto n_testcases = std::vector<testcase>{};
      for (auto h = 0; h < numHaplotypes; h++)
        n_testcases.push_back({.rslen = readLength,
                               .haplen = haplotypeLengths[h],
                               .q = readQuals,
                               .i = insGops,
                               .d = delGops,
                               .c = gapConts,
                               .hap = haplotypes[h],
                               .rs = reads});
      m_testcases.push_back(std::move(n_testcases));
    }
    return m_testcases;  
  }

  inline auto initNative(bool use_double = false, int max_threads = 64) {
    SPDLOG_DEBUG("Enter");

    g_use_double = use_double;
#ifdef _OPENMP
    const auto avail_threads = static_cast<int>(omp_get_max_threads());
    const auto req_threads = max_threads;
    g_max_threads = std::min(req_threads, avail_threads);

    SPDLOG_DEBUG("Available threads: {}", avail_threads);
    SPDLOG_DEBUG("Requested threads: {}", req_threads);
    if (req_threads > avail_threads) {
      SPDLOG_DEBUG("Using {} available threads, but {} were requested", g_max_threads,
        req_threads);
    } else {
      SPDLOG_DEBUG("Using {} threads", g_max_threads);
    }
#else
  if (max_threads != 1) {
    SPDLOG_DEBUG("Ignoring request for {} threads; not using OpenMP implementation",
      max_threads);
  }
#endif

    // enable FTZ
    if (_MM_GET_FLUSH_ZERO_MODE() != _MM_FLUSH_ZERO_ON) {
      SPDLOG_DEBUG("Flush-to-zero (FTZ) is enabled when running PairHMM");
    }
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

    g_compute_full_prob_float = &compute_full_prob_avxs;
    g_compute_full_prob_double = &compute_full_prob_avxd;

    // init convert char table
    ConvertChar::init();
    SPDLOG_DEBUG("Exit");
  }

  inline auto computeLikelihoodsNative(const std::vector<SamRecord<>>& readDataArray,
                           const std::vector<Haplotype>& haplotypeDataArray,
                           std::vector<std::vector<double>>& likelihoodArray) {
    SPDLOG_DEBUG("Enter");
                            
    //==================================================================
    // get data
    auto testcases = getData(readDataArray, haplotypeDataArray);
                            
    //==================================================================
    // calcutate pairHMM
                            
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(g_max_threads)
#endif
    for (auto i = 0; i < testcases.size(); i++) {
      for (auto j = 0; j < testcases[0].size(); j++) {
        auto result_final = 0.0;

        const auto result_float
          = g_use_double ? 0.0f : g_compute_full_prob_float(&testcases[i][j]);
        if (result_float < MIN_ACCEPTED) {
          const auto result_double = g_compute_full_prob_double(&testcases[i][j]);
          result_final
            = std::log10(result_double) - g_ctxd.LOG10_INITIAL_CONSTANT;
        } else {
          result_final = ::log10f(result_float) - g_ctxf.LOG10_INITIAL_CONSTANT;
        }
        likelihoodArray[i][j] = result_final;
        SPDLOG_DEBUG("result = {}", result_final);
      }
    }

    //==================================================================
    // release data
    SPDLOG_DEBUG("Exit");
  }

 public:
  auto
  compute_likelihoods(const std::vector<Haplotype>& haplotypeDataArray,
                      std::vector<SamRecord<>>& readDataArray) {
    initNative();
    auto likelihoodArray = std::vector(
      readDataArray.size(), std::vector(haplotypeDataArray.size(), 0.0));
    computeLikelihoodsNative(readDataArray, haplotypeDataArray,
                             likelihoodArray);
    PairHMM::normalize_likelihoods(likelihoodArray);
    PairHMM::filter_poorly_modeled_reads(readDataArray, likelihoodArray);
    return likelihoodArray;
  }
};

}  // namespace biovoltron
