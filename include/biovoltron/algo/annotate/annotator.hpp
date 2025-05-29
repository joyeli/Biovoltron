#pragma once

#include <biovoltron/algo/annotate/tree/interval_tree.hpp>
#include <unordered_map>

namespace biovoltron {

/**
 * @ingroup annotate
 * @brief TBA
 */
template<typename Data>
struct Annotator {
 private:
  std::unordered_map<std::string, IntervalTree<Data>> trees;

 public:
  template<std::convertible_to<Interval> L>
  auto
  insert_at(const Data& data, const L& location) {
    const auto interval = Interval{location};
    const auto key = std::string{interval.chrom + interval.strand};
    trees[key].insert(interval.begin, interval.end, data);
  }

  auto
  insert(const Data& data) requires std::convertible_to<Data, Interval> {
    const auto interval = Interval{data};
    const auto key = std::string{interval.chrom + interval.strand};
    trees[key].insert(interval.begin, interval.end, data);
  }

  auto
  index() {
    for (auto& [chrom, tree] : trees) { tree.index(); }
  }

  auto
  find(const Interval& interval) const {
    const auto key = std::string{interval.chrom + interval.strand};
    if (trees.contains(key)) {
      const auto& tree = trees.find(key)->second;
      return tree.find(interval.begin, interval.end);
    } else
      return std::vector<Data>{};
  }
};

}  // namespace biovoltron
