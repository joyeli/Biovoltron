#pragma once

#include <algorithm>
#include <cmath>
#include <ranges>

namespace biovoltron {

/**
 * @brief Extracts optimal and suboptimal alignment scores from a sorted score list.
 * @param scores List of scores (scores must be sorted in descending order).
 * @param diff Maximum allowed score difference for suboptimal count(default = 5).
 * @return Tuple {optimal score, suboptimal score, number of alignments with suboptimal score}.
 * @details
 *  - If only one score exists, suboptimal score and count are zero.
 *  - Stops counting suboptimal scores when score < (suboptimal - diff).
 */
constexpr auto
get_opt_subopt_count(auto scores, int diff = 5) noexcept
  -> std::tuple<int, int, int> {
  if (scores.empty())
    return {};
  if (scores.size() == 1)
    return {scores.front(), 0, 0};

  const auto opt_score = scores.front();
  const auto sub_score = scores[1];
  const auto min_score = sub_score - diff;

  const auto remains = scores | std::views::drop(1);
  const auto it = std::ranges::find_if(
    remains, [min_score](const auto score) { return score < min_score; });

  const auto sub_cnt = std::ranges::distance(remains.begin(), it);
  return {opt_score, sub_score, sub_cnt};
}

/**
 * @brief Container for alignment scoring metrics used in MAPQ computation.
 */
struct MemAln {
  int score{};      ///< Local best Smith–Waterman score.
  int score2{};     ///< Local second-best SW score(within same region).
  int sub_score{};  ///< Global second-best SW score(across all regions).
  int align_len{};  ///< Alignment length in bases.
  int sub_n{};      ///< Number of alignments with global second-best score.
  float frac_rep{}; ///< Fraction of repetitive k-mers in alignment region.
};

/**
 * @brief Approximates single-end mapping quality (MAPQ).
 * @param aln Alignment scoring metrics (see @ref MemAln).
 * @return Estimated MAPQ (0–60).
 * @details
 *  - Uses identity (fraction of matches) and read length scaling.
 *  - Penalizes multiple equally good hits (`sub_n`) and repetitive content (`frac_rep`).
 *  - Clamps output to range [0, 60].
 */
static auto
mem_approx_mapq_se(MemAln aln) {
  // Decompose alignment metrics
  const auto [score, csub, sub_score, l, sub_n, frac_rep] = aln;

  // Use sub_score if available, otherwise default to 20
  auto sub = sub_score ? sub_score : 20;
  sub = csub > sub ? csub : sub;
  // If suboptimal score is not less than best score, mapping is unreliable
  if (sub >= score)
    return 0;

  // Calculate identity and scaling factor
  const auto identity = 1. - (l - score) / 5. / l;
  // Scale by read length; short reads use 1, longer reads use log scaling
  auto tmp = l < 50 ? 1. : 3.912 / std::log(l);
  tmp *= identity * identity;

  // Main MAPQ formula: higher score difference and identity yield higher MAPQ
  auto mapq = static_cast<int>(6.02 * (score - sub) * tmp * tmp + .499);
  // Penalize for multiple suboptimal alignments
  if (sub_n > 0)
    mapq -= static_cast<int>(4.343 * std::log(sub_n + 1) + .499);

  // Clamp MAPQ to [0, 60]
  if (mapq > 60) mapq = 60;
  if (mapq < 0) mapq = 0;

  // Penalize for repetitive content
  mapq = static_cast<int>(mapq * (1. - frac_rep) + .499);
  return mapq;
}

/**
 * @brief Computes raw MAPQ from a score difference.
 * @param diff Difference between best and suboptimal alignment score.
 * @return MAPQ estimate.
 */
static auto
raw_mapq(int diff) noexcept {
  return static_cast<int>(6.02 * diff + .499);
}

/**
 * @brief Estimates paired-end mapping quality for both mates.
 * @param p0 Alignment metrics for read 1.
 * @param p1 Alignment metrics for read 2.
 * @param score_un Best unpaired alignment score.
 * @param o Best paired alignment score.
 * @param subo Suboptimal paired alignment score.
 * @param n_sub Number of suboptimal paired alignments.
 * @return Pair {MAPQ_read1, MAPQ_read2}, each clamped to [0, 60].
 * @details
 *  - Combines paired-end penalty scaling and single-end MAPQ.
 *  - Penalizes repetitive content and multiple suboptimal hits.
 *  - Adjusts single-end MAPQ if paired-end MAPQ is higher.
 */
static auto
mem_mapq_pe(MemAln p0, MemAln p1, int score_un, int o, int subo, int n_sub) {
  // Paired-end MAPQ calculation: combine score difference, repetitive content, and suboptimal hits
  subo = subo > score_un ? subo : score_un;
  auto q_pe = raw_mapq(o - subo);
  if (n_sub > 0)
    q_pe -= static_cast<int>(4.343 * std::log(n_sub + 1) + .499);
  if (q_pe < 0) q_pe = 0;
  if (q_pe > 60) q_pe = 60;

  // Penalize for repetitive content in both mates
  q_pe = static_cast<int>(q_pe * (1. - .5 * (p0.frac_rep + p1.frac_rep)) + .499);

  // Compute single-end MAPQ for each mate
  auto q_se0 = mem_approx_mapq_se(p0);
  auto q_se1 = mem_approx_mapq_se(p1);

  // Adjust single-end MAPQ if paired-end MAPQ is higher
  q_se0 = q_se0 > q_pe ? q_se0 : q_pe < q_se0 + 40 ? q_pe : q_se0 + 40;
  q_se1 = q_se1 > q_pe ? q_se1 : q_pe < q_se1 + 40 ? q_pe : q_se1 + 40;

  // Clamp single-end MAPQ to raw MAPQ based on score difference
  q_se0 = q_se0 < raw_mapq(p0.score - p0.score2) ? q_se0 : raw_mapq(p0.score - p0.score2);
  q_se1 = q_se1 < raw_mapq(p1.score - p1.score2) ? q_se1 : raw_mapq(p1.score - p1.score2);

  if (q_se0 > 60) q_se0 = 60;
  if (q_se1 > 60) q_se1 = 60;

  return std::pair{q_se0, q_se1};
}

}  // namespace biovoltron
