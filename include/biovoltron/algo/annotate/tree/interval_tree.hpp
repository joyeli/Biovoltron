#pragma once

#include <algorithm>
#include <biovoltron/utility/interval.hpp>
#include <stack>

namespace biovoltron {

/**
 * @ingroup annotate
 * @brief
 * Implement a  implicite augmented tree, implicite means the tree
 * is stored in a array. (based on cgrange)
 *
 * Augmented tree
 * A type of interval tree for retreving overlapping interval.
 *
 * Construction: O(NlogN)
 * Each node of the tree contains an interval and a integer, max,
 * recording the maximum end position among the subtree of current
 * node (including the current node).
 * Build a BST or self-balancing BST using begin position as key.
 * Update "max" in a buttom-up manner.
 *
 * Query: O(logN + M), where M is the number of hit.
 * Hit happens when
 *  1. query.begin < node.end
 *  2. query.end   > node.begin
 * Thus, we can eliminate
 *  1. the right subtree if node.begin > query.end
 *  2. both of the subtree if max < query.begin
 *
 * Implicite complete binary search tree
 * The tree is stored in an array. If we traverse the array, we will
 * be traversing the tree in in-order.
 *
 * Assume leaves are level 0, root are level K, then
 * 1. height: K+1
 * 2. total node: 2^(K+1) - 1
 * 3. at level k
 *  - k-th bit of the index is 0
 *  - lowest k-1 bits are all 1
 *  - first node at level k is indexed by 2^k - 1
 *  - a node with index x
 *    - left child: x - 2^(k-1)
 *    - right child: x + 2^(k-1)
 *    - is left child if k+1 bit is 0 (parent is x + 2^k), otherwise
 *      is right child (parent is x- 2^k)
 *    - there are 2^(k+1) nodes in the subtree (including current node)
 *    - left-most node of the subtree is x & ~(2^k - 1), i.e masking lower
 *      k bit to 0
 *
 * Imaginary node
 * If the tree is not complete, we imagin there are some imaginary node
 * to make the tree complete.
 */
template<typename Data>
struct IntervalTree {
 private:
  using index_type = std::uint32_t;

  struct Node {
    index_type begin{};
    index_type end{};
    Data data{};
    index_type max{};
  };

  struct Cell {
    index_type node{}, level{};
    bool left_processed{};
  };

  std::vector<Node> tree;

  index_type max_level;

  bool indexed;

  auto
  get_max_level() const {
    auto level = 1;
    for (; (1u << level) <= tree.size(); level++)
      ;
    return level - 1;
  }

  auto
  get_root_index() const {
    return (1u << get_max_level()) - 1;
  }

  auto
  get_subtree_range(index_type index, index_type level) {
    const auto leftmost_index = index >> level << level;
    const auto num_subtree_nodes = (1u << (level + 1)) - 1;
    const auto begin = tree.begin() + leftmost_index;
    const auto end = (begin + num_subtree_nodes) >= tree.end() ?
                       tree.end() :
                       begin + num_subtree_nodes;
    return std::tie(begin, end);
  }

  static auto
  left_child(index_type index, index_type level) {
    return index - (1u << (level - 1));
  }

  static auto
  right_child(index_type index, index_type level) {
    return index + (1u << (level - 1));
  }

  static auto
  is_left_child(index_type index, index_type level) {
    return ((index >> level) & 1) == 0;
  }

  static auto
  parent(index_type index, index_type level) {
    return is_left_child(index, level) ? index + (1u << level) :
                                         index - (1u << level);
  }

  auto
  update_max() {
    // set max for all leaves
    for (auto idx = 0; idx < tree.size(); idx += 2)
      tree[idx].max = tree[idx].end;

    // set max for imaginary node
    auto last_i = tree.size() - 1;
    auto last_max = tree[last_i].end;

    // set max level-by-level bottom-up
    auto level = 1;
    for (; (1u << level) <= tree.size(); level++) {
      const auto first_node = (1u << level) - 1;
      const auto step = (1u << (level + 1));
      for (auto idx = first_node; idx < tree.size(); idx += step) {
        const auto right_max = right_child(idx, level) < tree.size() ?
                                 tree[right_child(idx, level)].max :
                                 last_max;
        const auto left_max = tree[left_child(idx, level)].max;
        tree[idx].max = std::max({right_max, left_max, tree[idx].end});
      }
      // update max for imaginary node
      last_i = parent(last_i, level);
      if (last_i < tree.size() && tree[last_i].max > last_max)
        last_max = tree[last_i].max;
    }

    return level - 1;
  }

 public:
  auto
  insert(index_type begin, index_type end, const Data& data) {
    tree.emplace_back(Node{begin, end, data});
    indexed = false;
  }

  auto
  index() {
    if (!indexed) {
      std::ranges::sort(tree, {}, &Node::begin);
      max_level = update_max();
      indexed = true;
    }
  }

  auto
  find(index_type qbegin, index_type qend) const {
    if (!indexed)
      throw std::runtime_error("Try to find() before index()");

    auto results = std::vector<Data>{};
    auto jobs = std::stack<Cell>{};
    // push the root
    jobs.emplace((1u << max_level) - 1, max_level, false);
    while (!jobs.empty()) {
      auto job = jobs.top();
      jobs.pop();
      if (job.level <= 3) {
        // brute-force: traverse all nodes in subtree
        // TODO: better syntax?
        // auto [subtree_begin, subtree_end] = get_subtree_range(job.node,
        // job.level); std::transform(
        //   subtree_begin,
        //   subtree_end,
        //   std::back_inserter(results),
        //   [begin, end](const auto& node) { if (node.begin < end && node.end >
        //   begin) return node.data; });
        const auto begin = job.node >> job.level << job.level;
        const auto end = std::min((unsigned int)tree.size(), begin + (1u << (job.level + 1)) - 1);
        for (auto idx = begin; idx < end && tree[idx].begin < qend; idx++)
          if (qbegin < tree[idx].end)
            results.emplace_back(tree[idx].data);
      } else if (!job.left_processed) {
        job.left_processed = true;
        jobs.emplace(job);
        const auto left = left_child(job.node, job.level);
        if (left >= tree.size() || tree[left].max > qbegin)
          jobs.emplace(left, job.level - 1, false);
      } else if (job.node < tree.size() && tree[job.node].begin < qend) {
        if (qbegin < tree[job.node].end)
          results.emplace_back(tree[job.node].data);
        jobs.emplace(right_child(job.node, job.level), job.level - 1, false);
      }
    }
    return results;
  }
};

}  // namespace biovoltron
