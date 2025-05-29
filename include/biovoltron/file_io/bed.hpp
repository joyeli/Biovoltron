#pragma once

#include <biovoltron/file_io/core/header.hpp>
#include <biovoltron/file_io/core/record.hpp>
#include <biovoltron/utility/interval.hpp>

namespace biovoltron {

/**
 * @ingroup file_io
 * @brief TBA
 */
struct BedHeader : Header {
  constexpr static auto START_SYMBOLS = std::array{"browser", "track", "#"};
};

/**
 * @ingroup file_io
 * @brief TBA
 */
struct BedRecord : HeaderableRecord {
  BedHeader* header = nullptr;

  std::string chrom;
  std::uint32_t start{};
  std::uint32_t end{};
  std::string name;
  std::int32_t score{};
  char strand{};
  std::uint32_t thick_start{};
  std::uint32_t thick_end{};
  std::string item_rgb;
  std::uint32_t block_count{};
  std::string block_sizes;
  std::string block_starts;

  auto
  operator<=>(const BedRecord& other) const noexcept {
    return std::tie(chrom, start, end)
       <=> std::tie(other.chrom, other.start, other.end);
  }

  operator auto() const { return Interval{chrom, start, end, strand}; }
};

/**
 * @ingroup file_io
 * @brief TBA
 */
struct BedGraphRecord : HeaderableRecord {
  BedHeader* header = nullptr;

  std::string chrom;
  int start{};
  int end{};
  float score{};

  auto
  operator<=>(const BedRecord& other) const noexcept {
    return std::tie(chrom, start, end)
       <=> std::tie(other.chrom, other.start, other.end);
  }
};

}  // namespace biovoltron
