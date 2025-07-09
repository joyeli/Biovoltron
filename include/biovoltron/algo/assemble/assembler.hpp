#pragma once

#include <spdlog/spdlog.h>
#include <concepts>
#include <iostream>

#include <biovoltron/file_io/all.hpp>
#include <biovoltron/utility/haplotype/haplotype.hpp>
#include <biovoltron/algo/assemble/graph/adapter_graph.hpp>
#include <biovoltron/algo/assemble/graph/haplotype_graph.hpp>

namespace biovoltron {

class Assembler {};

class HaplotypeAssembler : Assembler {
private:
  template<std::derived_from<Record> R, class RefType>
    requires std::same_as<
      RefType, std::conditional_t<R::encoded, istring_view, std::string_view>
    >
  auto assemble(const std::vector<R> reads,
                RefType ref,
                std::size_t kmer_size) const -> std::vector<Haplotype> {
    
    static constexpr auto MIN_UNIQUE_KMERS_COUNT_TO_DISCARD = 4096;

    if (ref.size() < kmer_size) {
      return {};
    }

    auto graph = HaplotypeGraph(kmer_size);
    graph.set_ref(ref);
    for (const auto& read : reads) {
      graph.set_read(read);
    }

    graph.build();

    if (graph.unique_kmers_count() > MIN_UNIQUE_KMERS_COUNT_TO_DISCARD) {
      SPDLOG_DEBUG("Not using kmer size of {} in assembler because it contains too much unique kmers", kmer_size);
      return {};
    }

    if (graph.has_cycles()) {
      SPDLOG_DEBUG("Not using kmer size of {} in assembler because it contains a cycle", kmer_size);
      return {};
    }

    SPDLOG_DEBUG("Using kmer size of {} in assembler", kmer_size);
    return graph.find_paths();
  }

public:
  template<std::derived_from<Record> R, class RefType>
    requires std::same_as<
      RefType, std::conditional_t<R::encoded, istring_view, std::string_view>
    >
  auto assemble(const std::vector<R>& seqs, RefType ref) const {
    static constexpr auto INITIAL_KMER_SIZE = 25;
    static constexpr auto KMER_SIZE_ITERATION_INCREASE = 10;
    static constexpr auto MAX_ITERATIONS_TO_ATTEMPT = 6;

    auto iter = 1;
    auto kmer_size = INITIAL_KMER_SIZE;
    while (iter < MAX_ITERATIONS_TO_ATTEMPT) {
      auto haplotypes = assemble(seqs, ref, kmer_size);
      if (!haplotypes.empty()) {
        return haplotypes;
      }
      iter++;
      kmer_size += KMER_SIZE_ITERATION_INCREASE;
    }
    return std::vector<Haplotype>{};
  }
};

class AdapterAssembler : Assembler {

private:
  struct Parameter {
    /* the rate that consider the assembled adpater is adapter */
    const double LOW_COMPLEXITY_RATE = 0.7;
    const std::size_t INIT_KMER_SIZE = 10;
    const std::size_t INCREASE_KMER_SIZE = 5;
    const std::size_t MAX_KMER_SIZE = 35;

    /* @@@ non-sensitive mode parameters @@@ */
    /* the maximum iterations that try to assemble adapter */
    const std::size_t MAX_ITERATIONS = 3;

    /* the minimum number for occurence of the assembled adapter */
    const std::size_t MINIMUM_OCCURENCE = 10;

    /* the minimum rate of occurence of the assembled adapter */
    double PRUNE_FACTOR = 0.03;
    
    /* @@@ sensitive mode parameters @@@ */
    const std::size_t SEN_MAX_ITERATIONS = 5;
    const std::size_t SEN_MINIMUM_OCCURENCE = 0;
  } param;

private:

  template<class SeqType>
  auto assemble(const std::vector<SeqType>& seqs, 
                const std::size_t kmer_size,
                const std::size_t minimum_occurance) noexcept {
  
    auto graph = AdapterGraph<SeqType>(kmer_size,
                                       minimum_occurance);
    graph.build(seqs);
    auto adapters = graph.get_adapters();
    return adapters;
  }

  /**
   * @brief determine a sequence is low complexity or not. Low complexity is
   * defined as the assembled sequence have more than LOW_COMPLEXITY_RATE 
   * protion base are same.
   * @tparam SeqType Sequence type, must be std::string or istirng
   * @param seq the sequnece
   * @return `bool`
   */
  template<class SeqType>
  auto low_complexity(const SeqType& seq) const -> bool {
    auto cnt = std::array<int, 4>{};
    for (auto& base : seq) {
      if constexpr (std::same_as<SeqType, std::string>) {
        cnt[Codec::to_int(base)] += 1;
      } else {
        cnt[base] += 1;
      }
    }
    return std::ranges::any_of(cnt, [this, &seq](auto v) {
      return v >= seq.size() * param.LOW_COMPLEXITY_RATE;
    });
  };

public:

  /**
   * @brief Set the prune factor for adapter assembler
   *
   * @param new_factor
   * @return None
   */
  auto set_prune_factor(const double new_factor) noexcept {
    if (new_factor < 0 || 1 < new_factor) {
      SPDLOG_ERROR("The value of prune factor must inside [0 - 1], use default value {}", param.PRUNE_FACTOR);
    } else {
      param.PRUNE_FACTOR = new_factor;
    }
  }


  /**
   * @brief assemble adapter from given sequences
   * @tparam Seq
   * @param seqs sequence for assemble adapter
   * @param sensitive assembler adapter in sensitive mode
   * @return std::string if Seq is std::string, istring otherwise.
   */
  template<class Seq>
  auto assemble(const std::vector<Seq>& seqs, bool sensitive = false) noexcept {
    /**
     * default kmer size is set to 10 and the pruning factor, which is the
     * minimum percentage of k-mer occurence, is set to 0.03 by default.
     * We also add a constraint on the minimum number of occurence of a k-mer
     * to 10
     */
    const auto max_iter = 
      sensitive ? param.SEN_MAX_ITERATIONS : param.MAX_ITERATIONS;
    auto kmer_size = param.INIT_KMER_SIZE;
    const auto max_kmer_size = param.MAX_KMER_SIZE;

    auto iter = 1u;
    auto is_low_complexity = false;
    while (iter <= max_iter && kmer_size <= max_kmer_size) {

      /**
       * prune factor will decrease if cannot assemble proper adapter when
       * enter further iterations
       */
      const auto pruning_factor = param.PRUNE_FACTOR / iter;
      const auto minimum_occurence = std::max(
        static_cast<std::size_t>(std::ceil(seqs.size() * pruning_factor)),
        sensitive ? param.SEN_MINIMUM_OCCURENCE : param.MINIMUM_OCCURENCE
      );
      SPDLOG_DEBUG("Run iter {}", iter);
      SPDLOG_DEBUG("kmer size = {}", kmer_size);
      SPDLOG_DEBUG("Pruning factor = {}", pruning_factor);
      SPDLOG_DEBUG("Minimum occurance = {}", minimum_occurence);

      spdlog::info("Try to detect adapter with prune factor = {:.4f}", pruning_factor);
      auto adapters = assemble( 
        seqs,
        kmer_size,
        minimum_occurence
      );
      
      /* cannot assemble proper adapters */
      if (adapters.empty()) {
        spdlog::info("Failed :(");
        iter++;
        continue;
      }
      auto& adapter = adapters[0];
      /* the assembled adpater have low complexity */
      if (low_complexity(adapter)) {
        spdlog::info("The assembled adapter has low complexity");
        is_low_complexity = true;
        kmer_size += 5;
        continue;
      }
      /**
       * The output length of final adapter is set to 15 for a low-complexity
       * adapter since the bases near the end of the adapter are more likely
       * error-prune due to the nature of the de Bruijn graph
       */
      if (is_low_complexity) {
        adapter.resize(std::min(adapter.size(), 15ul));
      } else {
        adapter.resize(std::min(adapter.size(), 32ul));
      }
      spdlog::info("Success: Assembled adapter = {}", adapter);
      return adapter;
    }
    if constexpr (std::same_as<Seq, istring> || std::same_as<Seq, istring_view>) {
      return istring{};
    } else {
      return std::string{};
    }
  }
};

}  // namespace biovoltron
