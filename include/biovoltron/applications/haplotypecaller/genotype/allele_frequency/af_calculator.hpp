#pragma once

#include <biovoltron/math/math_utils.hpp>
#include <biovoltron/utility/genotype/genotype.hpp>
#include <biovoltron/utility/range/range_utils.hpp>
#include <ranges>

namespace biovoltron {

struct AFCalculator {
  constexpr static auto REF_PSEUDOCOUNT = 10.0;
  constexpr static auto SNP_PSEUDOCOUNT = 0.01;
  constexpr static auto INDEL_PSEUDOCOUNT = 0.00125;
  constexpr static auto THRESHOLD_FOR_ALLELE_COUNT_CONVERGENCE = 0.1;
  constexpr static auto HOM_REF_GENOTYPE_INDEX = 0;
  constexpr static auto SPAN_DEL = "*";

  static auto
  log10_normalized_genotype_posteriors(
    std::span<const double> log10_genotype_likelihoods,
    std::span<const double> log10_allele_frequencies,
    std::span<const Genotype> genotypes) {
    auto log10_posteriors = std::vector(genotypes.size(), 0.0);
    for (auto i = 0; i < log10_posteriors.size(); i++) {
      const auto [allele1, allele2] = genotypes[i];
      const auto log10_combination_count
        = allele1 == allele2 ? 0 : std::log10(2);
      const auto sum
        = log10_allele_frequencies[allele1] + log10_allele_frequencies[allele2];
      log10_posteriors[i]
        = log10_combination_count + log10_genotype_likelihoods[i] + sum;
    }
    return MathUtils::normalize_log10(log10_posteriors);
  }

  static auto
  effective_allele_counts(std::span<const double> log10_genotype_likelihoods,
                          std::span<const double> log10_allele_frequencies,
                          std::span<const Genotype> genotypes) {
    const auto log10_genotype_posteriors = log10_normalized_genotype_posteriors(
      log10_genotype_likelihoods, log10_allele_frequencies, genotypes);

    auto log10_result = std::vector(log10_allele_frequencies.size(),
                                    std::numeric_limits<double>::lowest());
    for (auto i = 0; i < genotypes.size(); i++) {
      const auto [allele1, allele2] = genotypes[i];
      log10_result[allele1] = MathUtils::log10_sum_log10(
        log10_result[allele1], log10_genotype_posteriors[i]);
      log10_result[allele2] = MathUtils::log10_sum_log10(
        log10_result[allele2], log10_genotype_posteriors[i]);
    }
    for (auto& x : log10_result) x = std::pow(10, x);
    return log10_result;
  }

  static auto
  genotype_indices_with_only_ref_and_span_del(
    std::span<const std::string> alleles, std::span<const Genotype> genotypes) {
    const auto span_del_index = RangeUtils::index_of(alleles, SPAN_DEL);
    auto non_variant_indices = std::vector{0};
    non_variant_indices.push_back(
      RangeUtils::index_of(genotypes, Genotype{0, span_del_index}));
    return non_variant_indices;
  }

  static auto
  calculate(std::span<const double> log10_genotype_likelihoods,
            std::span<const std::string> alleles,
            std::span<const Genotype> genotypes) {
    const auto num_alleles = alleles.size();
    auto prior_pseudocounts = std::vector(num_alleles, 0.0);
    prior_pseudocounts.front() = REF_PSEUDOCOUNT;
    for (auto i = 1; i < num_alleles; i++)
      prior_pseudocounts[i] = alleles[i].size() == alleles[0].size() ?
                                SNP_PSEUDOCOUNT :
                                INDEL_PSEUDOCOUNT;

    auto log10_allele_frequencies
      = std::vector(num_alleles, -std::log10(num_alleles));

    auto allele_counts = std::vector(num_alleles, 0.0);
    for (auto allele_counts_maximum_difference
         = std::numeric_limits<double>::max();
         allele_counts_maximum_difference
         > THRESHOLD_FOR_ALLELE_COUNT_CONVERGENCE;) {
      const auto new_allele_counts = effective_allele_counts(
        log10_genotype_likelihoods, log10_allele_frequencies, genotypes);

      const auto allele_counts_difference = RangeUtils::binary_transform(
        allele_counts, new_allele_counts, std::minus{});
      allele_counts_maximum_difference = std::ranges::max(
        allele_counts_difference
        | std::views::transform([](auto x) { return std::abs(x); }));

      allele_counts = new_allele_counts;
      const auto posterior_pseudocounts = RangeUtils::binary_transform(
        prior_pseudocounts, allele_counts, std::plus{});

      // first iteration uses flat prior in order to avoid local minimum where
      // the prior + no pseudocounts gives such a low effective allele frequency
      // that it overwhelms the genotype likelihood of a real variant basically,
      // we want a chance to get non-zero pseudocounts before using a prior
      // that's biased against a variant
      log10_allele_frequencies
        = MathUtils::dirichlet_log10_mean_weights(posterior_pseudocounts);
    }

    const auto log10_genotype_posteriors = log10_normalized_genotype_posteriors(
      log10_genotype_likelihoods, log10_allele_frequencies, genotypes);

    auto log10p_of_zero_counts_by_allele = std::vector(num_alleles, 0.0);
    auto log10p_no_variant = 0.0;
    auto log10_absent_posteriors
      = std::vector<std::vector<double>>(num_alleles);
    if (std::ranges::find(alleles, SPAN_DEL) == alleles.end())
      log10p_no_variant += log10_genotype_posteriors[HOM_REF_GENOTYPE_INDEX];
    else {
      const auto non_variant_indices
        = genotype_indices_with_only_ref_and_span_del(alleles, genotypes);

      auto non_variant_log10_posteriors = std::vector<double>{};
      for (const auto n : non_variant_indices)
        non_variant_log10_posteriors.push_back(log10_genotype_posteriors[n]);

      log10p_no_variant += std::min(
        0.0, MathUtils::log10_sum_log10(non_variant_log10_posteriors));
    }

    if (num_alleles != 2) {
      for (auto i = 0; i < genotypes.size(); i++) {
        const auto log10_genotype_posterior = log10_genotype_posteriors[i];
        const auto [allele1, allele2] = genotypes[i];
        for (auto n = 0; n < log10_absent_posteriors.size(); n++)
          if (n != allele1 && n != allele2)
            log10_absent_posteriors[n].push_back(log10_genotype_posterior);
      }

      auto log10p_no_allele = std::vector<double>{};
      for (const auto& buffer : log10_absent_posteriors)
        log10p_no_allele.push_back(
          std::min(0.0, MathUtils::log10_sum_log10(buffer)));

      for (auto i = 0; i < log10p_of_zero_counts_by_allele.size(); i++)
        log10p_of_zero_counts_by_allele[i] += log10p_no_allele[i];
    }

    if (num_alleles == 2)
      log10p_of_zero_counts_by_allele[1] = log10p_no_variant;

    auto integer_allele_counts = std::vector<int>{};
    for (const auto x : allele_counts)
      integer_allele_counts.push_back(std::round(x));

    auto mle_counts_and_log10p_ref_by_allele
      = std::vector<std::pair<std::string, std::pair<int, double>>>{};
    for (auto i = 1; i < num_alleles; i++)
      mle_counts_and_log10p_ref_by_allele.emplace_back(
        alleles[i], std::pair{integer_allele_counts[i],
                              log10p_of_zero_counts_by_allele[i]});

    return std::pair{log10p_no_variant,
                     std::move(mle_counts_and_log10p_ref_by_allele)};
  }

  static auto
  calculate_single_sample_biallelic_non_ref_posterior(
    std::array<double, 3> log10_genotype_likelihoods) {
    if (std::ranges::max_element(log10_genotype_likelihoods)
        == log10_genotype_likelihoods.cbegin())
      return 0.0;
    auto log10_unnormalized_posteriors = std::array<double, 3>{};
    for (auto i = 0; i < log10_unnormalized_posteriors.size(); i++) {
      log10_unnormalized_posteriors[i]
        = log10_genotype_likelihoods[i]
          + MathUtils::log10_binomial_coefficient(2, i)
          + MathUtils::log_to_log10(std::lgamma(i + SNP_PSEUDOCOUNT)
                                    + std::lgamma(2 - i + REF_PSEUDOCOUNT));
    }

    if (std::ranges::max_element(log10_unnormalized_posteriors)
        == log10_unnormalized_posteriors.cbegin())
      return 0.0;
    return 1
           - MathUtils::normalize_from_log10_to_linear_space(
             log10_unnormalized_posteriors)[0];
  }
};

}  // namespace biovoltron
