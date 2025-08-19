#pragma once
#include <iostream>
#include <biovoltron/file_io/sam.hpp>

namespace biovoltron {

  /**
 * @ingroup utility
 * @brief Filter reads with low mapping quality.
 */
struct MappingQualityReadFilter {
  /**
   * Minimum acceptable mapping quality score for filtering.
   */
  static constexpr auto MIN_MAPPING_QUALITY_SCORE = 20;

  /**
   * Filters out reads with low mapping quality.
   *
   * @tparam Encoded Whether the read is encoded (default is false).
   * @param record The SAM record to check.
   * @return True if the read should be filtered out, false otherwise.
   */
  template<bool Encoded = false>
  auto
  operator()(const SamRecord<Encoded>& record) const noexcept {
    return record.mapq < MIN_MAPPING_QUALITY_SCORE;
  }
};

/**
 * @ingroup utility
 * @brief Filter duplicate reads.
 */
struct DuplicateReadFilter {
  /**
   * Filters out duplicate reads.
   * 
   * @tparam Encoded Whether the read is encoded (default is false).
   * @param record The SAM record to check.
   * @return True if the read is a duplicate, false otherwise.
   */
  template<bool Encoded = false>
  auto
  operator()(const SamRecord<Encoded>& record) const noexcept {
    return record.duplicate_read();
  }
};

/**
 * @ingroup utility
 * @brief Filter secondary alignments.
 */
struct SecondaryAlignmentReadFilter {
  /**
   * Filters out secondary alignments.
   *
   * @tparam Encoded Whether the read is encoded (default is false).
   * @param record The SAM record to check.
   * @return True if the read is a secondary alignment, false otherwise.
   */
  template<bool Encoded = false>
  auto
  operator()(const SamRecord<Encoded>& record) const noexcept {
    return record.secondary_alignment();
  }
};

/**
 * @ingroup utility
 * @brief Filter reads shorter than the minimum length.
 */
struct MinimumLengthReadFilter {
  /**
   * Minimum read length after trimming.
   */
  static constexpr auto MINIMUM_READ_LENGTH_AFTER_TRIMMING = 10;

  /**
   * Filters out reads shorter than the minimum length.
   *
   * @tparam Encoded Whether the read is encoded (default is false).
   * @param record The SAM record to check.
   * @return True if the read is shorter than the minimum length, false otherwise.
   */
  template<bool Encoded = false>
  auto
  operator()(const SamRecord<Encoded>& record) const noexcept {
    return record.size() < MINIMUM_READ_LENGTH_AFTER_TRIMMING;
  }
};

/**
 * @ingroup utility
 * @brief Filter reads whose mate is on a different contig.
 */
struct MateOnSameContigReadFilter {
  /**
   * @brief Filters out reads whose mate is on a different contig.
   *
   * @tparam Encoded Whether the read is encoded (default is false).
   * @param record The SAM record to check.
   * @return True if the read's mate is on a different contig, false otherwise.
   */
  template<bool Encoded = false>
  auto
  operator()(const SamRecord<Encoded>& record) const noexcept {
    return record.rnext != "=";
  }
};

}  // namespace biovoltron