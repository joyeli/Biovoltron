#pragma once

#include <array>
#include <cmath>

namespace biovoltron {

/**
 * @ingroup utility
 * @brief Utilities for quality score conversions.
 */
struct QualityUtils {

  /**
   * @brief ASCII offset for FASTQ quality ('!').
   */
  static constexpr char ASCII_OFFSET = '!';

  /**
   * @brief Convert a quality score character to its error probability.
   *
   * @param qual Quality score character (e.g., '!' for Q = 0).
   * @return Error probability corresponding to the quality score.
   */
  static double
  qual_to_error_prob(char qual) {
    return qual_to_error_prob_cache[qual];
  }

  /**
   * @brief Convert a quality score to its error probability in log10 scale.
   * 
   * @param qual Quality score.
   * @return Error probability in log10 scale corresponding to the quality score.
   */

  static auto
  qual_to_error_prob_log10(double qual) {
    return qual / -10.0;
  }

  /**
   * @brief Convert a quality score to its probability in log10 scale.
   *
   * @param qual Quality score.
   * @return Probability in log10 scale corresponding to the quality score.
   */
  static auto
  qual_to_prob_log10(double qual) {
    return std::log10(1 - qual_to_error_prob(qual));
  }

  /**
   * @brief Convert an error rate to a quality score.
   *
   * @param error_rate Error rate (0.0 <= error_rate <= 1.0).
   * @return Quality score corresponding to the error rate.
   */
  static auto
  phred_scale_error_rate(double error_rate) {
    return -10.0 * std::log10(error_rate);
  }

/**
 * @brief Lookup table for quality scores to error probabilities.
 *
 * This array stores `10^(-Q/10)` for Q = 0..127,
 * where Q is the ASCII value minus `ASCII_OFFSET`.
 */
 private:
  constexpr static auto qual_to_error_prob_cache = [] {
    auto qual_to_error_prob_cache = std::array<double, 128>{};
    for (auto qual = 0; qual < qual_to_error_prob_cache.size(); qual++)
      qual_to_error_prob_cache[qual] = std::pow(10, qual / -10.0);
    return qual_to_error_prob_cache;
  }();
};

}  // namespace biovoltron