#pragma once

#include <algorithm>
#include <span>

namespace biovoltron {

class DepthPerAllele {
  constexpr static auto LOG_10_INFORMATIVE_THRESHOLD = 0.2;

  static auto
  get_informative_alleles(std::span<const std::vector<double>> likelihoods) {
    auto result = std::vector<int>{};
    for (auto alleles : likelihoods) {
      const auto i = std::ranges::max_element(alleles) - alleles.begin();
      std::ranges::sort(alleles, std::ranges::greater{});
      const auto confidence = alleles[0] - alleles[1];
      if (confidence > LOG_10_INFORMATIVE_THRESHOLD)
        result.push_back(i);
    }
    return result;
  }

 public:
  static auto
  annotate(std::span<const std::vector<double>> likelihoods) {
    auto counts = std::vector(likelihoods.front().size(), 0);
    for (const auto allele : get_informative_alleles(likelihoods))
      counts[allele]++;
    return counts;
  }
};

}  // namespace biovoltron