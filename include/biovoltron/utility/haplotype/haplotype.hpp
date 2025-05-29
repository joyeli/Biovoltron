#pragma once

#include <biovoltron/file_io/cigar.hpp>
#include <biovoltron/utility/interval.hpp>
#include <biovoltron/utility/variant/variant.hpp>
#include <cmath>
#include <map>

namespace biovoltron {

struct Haplotype {
  std::string seq;
  Interval location;
  std::map<std::uint32_t, Variant> event_map;
  Cigar cigar;
  std::uint32_t align_begin_wrt_ref = 0;
  double score = std::numeric_limits<double>::lowest();
  int rank;

  auto
  size() const noexcept {
    return seq.size();
  }

  auto
  get_overlapping_events(std::size_t begin) const {
    auto events =
      std::ranges::subrange(event_map.begin(), event_map.upper_bound(begin))
      | std::views::values | std::views::filter([begin](const auto& var) {
          return var.location.end > begin;
        });
    return std::vector(events.begin(), events.end());
  }
};

}  // namespace biovoltron