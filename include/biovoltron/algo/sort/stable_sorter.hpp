#pragma once

#include <biovoltron/utility/istring.hpp>
#include <biovoltron/algo/sort/core/sorter.hpp>
#include <execution>

namespace biovoltron {

/**
 * @ingroup sort
 * @brief Stable suffix array construction using parallel execution.
 * @tparam size_type Integer type used for suffix array indices (default = std::uint32_t).
 * @details
 *  - Builds a suffix array (SA) where each element is the starting position
 *    of a suffix in the reference string, sorted in lexicographical order.
 *  - Suffix arrays are widely used for fast substring search and FM-index construction.
 *  - Uses `std::stable_sort` to preserve the original order of equal keys, which can be important
 *    for reproducibility and downstream algorithms that rely on stable ordering.
 *  - Parallelized with `std::execution::par_unseq` for performance on large inputs.
 */

template<typename size_type = std::uint32_t>

struct StableSorter {
  /// Alias for suffix array container type.
  using SA_t = std::vector<size_type>;

  /**
   * @brief Constructs a stable suffix array from the given reference string.
   * @param ref The reference string (encoded as @ref istring_view).
   * @param sort_len Optional limit on the length of substring comparisons.
   *                 If `istring_view::npos`, compares until the end of each suffix.
   *                 Shorter sort_len can speed up sorting but may produce non-unique order for equal prefixes.
   * @return Suffix array (vector of starting positions, including sentinel position ref.size()).
   * @details
   *  - Generates an initial index list `[0, 1, ..., N]`, where N = ref.size().
   *  - Sorts indices based on lexicographical order of the substrings starting at those positions.
   *  - Stable sort ensures relative order of equal substrings is preserved.
   *  - The resulting suffix array can be used for substring search and FM-index construction.
   */
  static auto
  get_sa(istring_view ref, std::size_t sort_len = istring_view::npos) {
    auto sa = std::vector<size_type>(ref.size() + 1);
    std::iota(sa.begin(), sa.end(), 0);
    std::stable_sort(std::execution::par_unseq, sa.begin(), sa.end(),
                     [ref, sort_len](auto i, auto j) {
                       return ref.substr(i, sort_len) < ref.substr(j, sort_len);
                     });
    return sa;
  }
};

}  // namespace biovoltron
