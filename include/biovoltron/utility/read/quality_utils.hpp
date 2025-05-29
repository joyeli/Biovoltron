#pragma once

#include <array>
#include <cmath>

namespace biovoltron {

struct QualityUtils {
  static constexpr char ASCII_OFFSET = '!';

  static double
  qual_to_error_prob(char qual) {
    return qual_to_error_prob_cache[qual];
  }

  static auto
  qual_to_error_prob_log10(double qual) {
    return qual / -10.0;
  }

  static auto
  qual_to_prob_log10(double qual) {
    return std::log10(1 - qual_to_error_prob(qual));
  }

  static auto
  phred_scale_error_rate(double error_rate) {
    return -10.0 * std::log10(error_rate);
  }

 private:
  constexpr static auto qual_to_error_prob_cache = [] {
    auto qual_to_error_prob_cache = std::array<double, 128>{};
    for (auto qual = 0; qual < qual_to_error_prob_cache.size(); qual++)
      qual_to_error_prob_cache[qual] = std::pow(10, qual / -10.0);
    return qual_to_error_prob_cache;
  }();
};

}  // namespace biovoltron