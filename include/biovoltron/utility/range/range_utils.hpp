#pragma once

#include <algorithm>
#include <cassert>
#include <ostream>
#include <ranges>
#include <span>

namespace biovoltron {

struct RangeUtils {
  static auto
  binary_transform(std::span<const double> a, std::span<const double> b,
                   auto op) {
    assert(a.size() == b.size());
    auto result = std::vector(a.begin(), a.end());
    std::ranges::transform(result, b, result.begin(), op);
    return result;
  }

  static auto
  index_of(const auto& r, const auto& value) {
    return std::ranges::distance(r.begin(), std::ranges::find(r, value));
  };

  static auto
  format_print(const std::ranges::range auto& r, std::ostream& os,
               std::string delim = ",") {
    for (auto d = std::string{}; const auto& elem : r) {
      os << d << elem;
      d = delim;
    }
  }

  template<std::default_initializable Comp, std::ranges::random_access_range R>
    requires std::sortable<std::ranges::iterator_t<R>, Comp>
  static auto
  second(R r) {
    assert(r.size() >= 2);
    std::ranges::sort(r, Comp{});
    return r[1];
  }
};

}  // namespace biovoltron