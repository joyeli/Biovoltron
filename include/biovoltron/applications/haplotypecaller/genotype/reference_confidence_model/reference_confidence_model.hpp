#pragma once

#include <biovoltron/math/math_utils.hpp>
#include <biovoltron/utility/istring.hpp>
#include <biovoltron/utility/read/quality_utils.hpp>

namespace biovoltron {

struct ReferenceConfidenceModel {
  constexpr static auto MIN_BASE_QUAL = 10;
  constexpr static auto LOG10_ONE_THIRD = -0.47712125472;
  constexpr static auto REF_MODEL_DELETION_QUAL = 30;
  constexpr static auto LOG10_PLOIDY = 0.30103;  // log10(2)

  static auto
  apply_pileup_element_ref_vs_non_ref_likelihood_and_count(
    ichar ref_base, std::array<double, 3>& genotype_likelihoods,
    ichar read_base, ichar qual) {
    auto reference_likelihood = 0.0, non_ref_likelihood = 0.0;
    // is_alt
    if (read_base != ref_base) {
      non_ref_likelihood = QualityUtils::qual_to_prob_log10(qual);
      reference_likelihood
        = QualityUtils::qual_to_error_prob_log10(qual) + LOG10_ONE_THIRD;
    } else {
      reference_likelihood = QualityUtils::qual_to_prob_log10(qual);
      non_ref_likelihood
        = QualityUtils::qual_to_error_prob_log10(qual) + LOG10_ONE_THIRD;
    }
    genotype_likelihoods[0] += reference_likelihood + LOG10_PLOIDY;
    genotype_likelihoods[2] += non_ref_likelihood + LOG10_PLOIDY;
    genotype_likelihoods[1] += MathUtils::approximate_log10_sum_log10(
      reference_likelihood, non_ref_likelihood);
    return genotype_likelihoods;
  }

 private:
  static auto
  calc_genotype_likelihoods_of_ref_vs_any(const std::vector<ichar>& read_pileup,
                                          const std::vector<ichar>& qual_pileup,
                                          ichar ref_base) {
    assert(std::ranges::all_of(qual_pileup,
                               [](auto qual) { return qual > MIN_BASE_QUAL; }));

    auto genotype_likelihoods = std::array<double, 3>{};
    for (auto i = 0; i < read_pileup.size(); i++)
      apply_pileup_element_ref_vs_non_ref_likelihood_and_count(
        ref_base, genotype_likelihoods, read_pileup[i], qual_pileup[i]);

    const auto denominator = read_pileup.size() * LOG10_PLOIDY;
    for (auto& value : genotype_likelihoods) value -= denominator;
    return genotype_likelihoods;
  }
};

}  // namespace biovoltron