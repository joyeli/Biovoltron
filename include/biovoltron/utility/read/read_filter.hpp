#pragma once

#include <biovoltron/file_io/sam.hpp>

namespace biovoltron {

struct MappingQualityReadFilter {
  static constexpr auto MIN_MAPPING_QUALITY_SCORE = 20;
  template<bool Encoded = false>
  auto
  operator()(const SamRecord<Encoded>& record) const noexcept {
    return record.mapq < MIN_MAPPING_QUALITY_SCORE;
  }
};

struct DuplicateReadFilter {
  template<bool Encoded = false>
  auto
  operator()(const SamRecord<Encoded>& record) const noexcept {
    return record.duplicate_read();
  }
};

struct SecondaryAlignmentReadFilter {
  template<bool Encoded = false>
  auto
  operator()(const SamRecord<Encoded>& record) const noexcept {
    return record.secondary_alignment();
  }
};

struct MinimumLengthReadFilter {
  static constexpr auto MINIMUM_READ_LENGTH_AFTER_TRIMMING = 10;
  template<bool Encoded = false>
  auto
  operator()(const SamRecord<Encoded>& record) const noexcept {
    return record.size() < MINIMUM_READ_LENGTH_AFTER_TRIMMING;
  }
};

struct MateOnSameContigReadFilter {
  template<bool Encoded = false>
  auto
  operator()(const SamRecord<Encoded>& record) const noexcept {
    return record.rnext != "=";
  }
};

}  // namespace biovoltron