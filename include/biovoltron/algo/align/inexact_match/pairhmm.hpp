#pragma once

#include <biovoltron/file_io/sam.hpp>
#include <biovoltron/utility/haplotype/haplotype.hpp>
#include <biovoltron/utility/read/quality_utils.hpp>

namespace biovoltron {

/**
 * @ingroup align
 *
 * @brief TBA
 */
struct PairHMM {
  enum {
    M_TO_M,
    M_TO_I,
    M_TO_D,
    I_TO_M,
    I_TO_I,
    D_TO_M,
    D_TO_D,
    TRANS_PROB_ARRAY_LENGTH
  };

  using TrasMatrx = std::array<double, TRANS_PROB_ARRAY_LENGTH>;
  static constexpr auto TRISTATE_CORRECTION = 3.0;
  static constexpr auto MAXIMUM_BEST_ALT_LIKELIHOOD_DIFFERENCE = -4.5;
  static constexpr auto EXPECTED_ERROR_RATE_PER_BASE = 0.02;
  static constexpr auto LOG10_QUALITY_PER_BASE = -4.0;
  static constexpr auto MAXIMUM_EXPECTED_ERROR_PER_READ = 2.0;
  static constexpr auto ORIGINAL_DEFAULT
    = TrasMatrx{0.9998, 0.0001, 0.0001, 0.9, 0.1, 0.9, 0.1};
  static inline const auto INITIAL_CONDITION = std::pow(2, 1020);
  static inline const auto INITIAL_CONDITION_LOG10
    = std::log10(INITIAL_CONDITION);

 private:
  auto
  initialize_priors(SamRecord<>& read, std::string_view haplotype,
                    auto& prior) const {
    for (auto i = 0; i < read.size(); i++) {
      const auto qual = read.qual[i];
      const auto x = read.seq[i];
      for (auto j = 0; j < haplotype.size(); j++) {
        const auto y = haplotype[j];
        prior[i + 1][j + 1]
          = (x == y || x == 'N' || y == 'N' ?
               1 - QualityUtils::qual_to_error_prob(qual) :
               (QualityUtils::qual_to_error_prob(qual) / TRISTATE_CORRECTION));
      }
    }
  }

  auto
  sub_compute_likelihood(SamRecord<>& read, std::string_view haplotype,
                         const TrasMatrx& trans, auto& M, auto& I, auto& D,
                         auto& prior) const {
    auto initial_value = INITIAL_CONDITION / haplotype.size();
    for (auto j = 0; j <= haplotype.size(); j++) D[0][j] = initial_value;

    initialize_priors(read, haplotype, prior);

    for (auto i = 1; i <= read.size(); i++) {
      for (auto j = 1; j <= haplotype.size(); j++) {
        M[i][j]
          = prior[i][j]
            * (M[i - 1][j - 1] * trans[M_TO_M] + I[i - 1][j - 1] * trans[I_TO_M]
               + D[i - 1][j - 1] * trans[D_TO_M]);
        I[i][j] = M[i - 1][j] * trans[M_TO_I] + I[i - 1][j] * trans[I_TO_I];
        D[i][j] = M[i][j - 1] * trans[M_TO_D] + D[i][j - 1] * trans[D_TO_D];
      }
    }

    auto final_sum_prob = 0.0;
    const auto end_i = read.size();
    for (auto j = 1; j <= haplotype.size(); j++)
      final_sum_prob += M[end_i][j] + I[end_i][j];
    return std::log10(final_sum_prob) - INITIAL_CONDITION_LOG10;
  }

 public:
  static auto
  normalize_likelihoods(std::vector<std::vector<double>>& log_likelihoods) {
    for (auto& likelihoods : log_likelihoods) {
        auto best_likelihood = std::ranges::max(likelihoods);
        auto cap_likelihood
          = best_likelihood + MAXIMUM_BEST_ALT_LIKELIHOOD_DIFFERENCE;
        for (auto& likelihood : likelihoods)
          if (likelihood < cap_likelihood)
            likelihood = cap_likelihood;
      }
  }

  static auto
  filter_poorly_modeled_reads(
    std::vector<SamRecord<>>& reads,
    std::vector<std::vector<double>>& log_likelihoods) {
    auto filtered_reads = std::vector<SamRecord<>>{};
    auto filtered_likelihoods = std::vector<std::vector<double>>{};
    for (auto i = 0; i < log_likelihoods.size(); i++) {
      auto best_likelihood = std::ranges::max(log_likelihoods[i]);
      auto likelihood_threshold
        = std::min(MAXIMUM_EXPECTED_ERROR_PER_READ,
                   std::ceil(reads[i].size() * EXPECTED_ERROR_RATE_PER_BASE))
          * LOG10_QUALITY_PER_BASE;
      if (best_likelihood < likelihood_threshold) {
        filtered_reads.push_back(std::move(reads[i]));
        filtered_likelihoods.push_back(std::move(log_likelihoods[i]));
      }
    }
    reads.swap(filtered_reads);
    log_likelihoods.swap(filtered_likelihoods);
  }

  auto
  compute_likelihoods(const std::vector<Haplotype>& haplotypes,
                      std::vector<SamRecord<>>& reads,
                      const TrasMatrx& trans = ORIGINAL_DEFAULT) const {
    auto get_padded_size = [](const auto& r) {
      return 1 + std::ranges::max_element(r, {}, [](const auto& elem) {
                   return elem.size();
                 })->size();
    };

    auto M = std::vector(get_padded_size(reads),
                         std::vector(get_padded_size(haplotypes), 0.0));
    auto I = M;
    auto D = M;
    auto prior = M;

    auto log_likelihoods
      = std::vector(reads.size(), std::vector(haplotypes.size(), 0.0));
    for (auto i = 0; i < reads.size(); i++)
      for (auto j = 0; j < haplotypes.size(); j++)
        log_likelihoods[i][j] = sub_compute_likelihood(
          reads[i], haplotypes[j].seq, trans, M, I, D, prior);
    normalize_likelihoods(log_likelihoods);
    filter_poorly_modeled_reads(reads, log_likelihoods);
    return log_likelihoods;
  }
};

}  // namespace biovoltron
