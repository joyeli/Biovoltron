#pragma once

#include <biovoltron/math/math_utils.hpp>
#include <biovoltron/utility/genotype/genotype_utils.hpp>
#include <biovoltron/utility/range/range_utils.hpp>

namespace biovoltron {

class AlleleUtils {
  static auto
  subsetted_pl_indices(std::span<const std::string> original_alleles,
                       std::span<const std::string> new_alleles,
                       std::span<const Genotype> genotypes) {
    auto new_allele_indices = std::vector<int>{};
    for (const auto& new_allele : new_alleles)
      new_allele_indices.push_back(
        RangeUtils::index_of(original_alleles, new_allele));

    const auto new_genotypes
      = GenotypeUtils::get_vcf_genotypes(new_allele_indices);

    auto result = std::vector<int>{};
    for (const auto& new_genotype : new_genotypes)
      result.push_back(RangeUtils::index_of(genotypes, new_genotype));
    return result;
  }

 public:
  static auto
  subset_alleles(std::span<const double> log10_genotype_likelihoods,
                 std::span<const std::string> original_alleles,
                 std::span<const std::string> alleles_to_keep,
                 std::span<const Genotype> genotypes) {
    const auto subsetted_likelihood_indices
      = subsetted_pl_indices(original_alleles, alleles_to_keep, genotypes);

    auto new_likelihoods = std::vector<double>{};
    for (const auto idx : subsetted_likelihood_indices)
      new_likelihoods.push_back(log10_genotype_likelihoods[idx]);
    new_likelihoods =

      MathUtils::scale_log_space_array_for_numerical_stability(new_likelihoods);

    return new_likelihoods;
  }
};

}  // namespace biovoltron