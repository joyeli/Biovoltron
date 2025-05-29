#pragma once

#include <biovoltron/utility/genotype/genotype.hpp>
#include <biovoltron/utility/range/range_utils.hpp>
#include <cmath>
#include <ranges>

namespace biovoltron {

struct GenotypeUtils {
  constexpr static auto MAX_ALLELE_COUNT = 7;

 private:
  static auto
  generate_vcf_genotypes(const std::ranges::range auto& alleles)
    -> std::vector<Genotype> {
    auto genotypes = std::vector<Genotype>{};
    for (const auto allele1 : alleles) {
      for (const auto allele2 : alleles) {
        if (allele2 > allele1)
          break;
        genotypes.emplace_back(allele2, allele1);
      }
    }
    return genotypes;
  }

  const static inline auto vcf_genotypes = [] {
    auto vcf_genotypes = std::vector<std::vector<Genotype>>{};
    for (auto i = 0; i <= MAX_ALLELE_COUNT; i++)
      vcf_genotypes.push_back(generate_vcf_genotypes(std::views::iota(0, i)));
    return vcf_genotypes;
  }();

  static auto
  generate_raw_genotypes(int num_alleles) -> std::vector<Genotype> {
    auto genotypes = std::vector<Genotype>{};
    for (auto allele1 = 0; allele1 < num_alleles; allele1++)
      for (auto allele2 = allele1; allele2 < num_alleles; allele2++)
        genotypes.emplace_back(allele1, allele2);
    return genotypes;
  }

  const static inline auto raw_genotypes = [] {
    auto raw_genotypes = std::vector<std::vector<Genotype>>{};
    for (auto i = 0; i <= MAX_ALLELE_COUNT; i++)
      raw_genotypes.push_back(generate_raw_genotypes(i));
    return raw_genotypes;
  }();

  const static inline auto raw_to_vcf_tables = [] {
    auto raw_to_vcf_tables = std::vector<std::vector<int>>{};
    for (auto i = 0; i < raw_genotypes.size(); i++) {
      auto& back = raw_to_vcf_tables.emplace_back();
      for (const auto& raw_genotype : raw_genotypes[i])
        back.push_back(RangeUtils::index_of(vcf_genotypes[i], raw_genotype));
    }
    return raw_to_vcf_tables;
  }();

 public:
  static auto
  get_vcf_genotypes(int num_alleles) {
    return std::span{vcf_genotypes.at(num_alleles)};
  }

  static auto
  get_vcf_genotypes(const std::ranges::range auto& alleles) {
    return generate_vcf_genotypes(alleles);
  }

  static auto
  get_raw_genotypes(int num_alleles) {
    return std::span{raw_genotypes.at(num_alleles)};
  }

  static auto
  get_genotype_size(int num_alleles) {
    return vcf_genotypes.at(num_alleles).size();
  }

  static auto
  get_allele_size(int num_genotypes) {
    return RangeUtils::index_of(
      vcf_genotypes
        | std::views::transform([](const auto& elem) { return elem.size(); }),
      num_genotypes);
  }

  static auto
  to_vcf_order(std::span<const double> raw_pls) {
    const auto allele_size = get_allele_size(raw_pls.size());
    const auto& table = raw_to_vcf_tables[allele_size];
    auto vcf_pls = std::vector(raw_pls.size(), 0.0);
    for (auto i = 0; i < raw_pls.size(); i++) vcf_pls[table[i]] = raw_pls[i];
    return vcf_pls;
  }

  static auto
  gls_to_pls(const std::vector<double>& gls) {
    auto pls = std::vector(gls.size(), 0);
    const auto adjust = std::ranges::max(gls);
    for (auto i = 0; i < pls.size(); i++)
      pls[i] = std::round(-10 * (gls[i] - adjust));
    return pls;
  }
};

}  // namespace biovoltron
