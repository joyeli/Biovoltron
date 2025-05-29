#pragma once

#include <biovoltron/file_io/cigar.hpp>
#include <cassert>
#include <limits>

namespace biovoltron {

/**
 * @ingroup align
 *
 * @brief TBA
 */
struct SmithWaterman {
  struct Parameters {
    int w_match;
    int w_mismatch;
    int w_open;
    int w_extend;
  };

  // match=1, mismatch = -1/3, gap=-(1+k/3)
  static constexpr auto ORIGINAL_DEFAULT = Parameters{3, -1, -4, -3};
  static constexpr auto STANDARD_NGS = Parameters{25, -50, -110, -6};
  static constexpr auto NEW_SW_PARAMETERS = Parameters{200, -150, -260, -11};
  static constexpr auto ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS
    = Parameters{10, -15, -30, -5};
  static constexpr auto MAX_MISMATCHES = 2;

  static auto
  well_match(std::string_view ref, std::string_view alt) {
    auto mismatch = 0;
    for (auto i = 0; mismatch <= MAX_MISMATCHES && i < ref.size(); i++)
      if (alt[i] != ref[i])
        mismatch++;
    return mismatch <= MAX_MISMATCHES;
  }

 private:
  static auto
  calculate_matrix(std::string_view ref, std::string_view alt,
                   std::vector<std::vector<int>>& score,
                   std::vector<std::vector<int>>& trace, Parameters params) {
    const auto row_size = score.size();
    const auto col_size = score[0].size();

    auto gap_size_down = std::vector(col_size + 1, 0);
    auto best_gap_down
      = std::vector(col_size + 1, std::numeric_limits<int>::min() / 2);
    auto gap_size_right = std::vector(row_size + 1, 0);
    auto best_gap_right
      = std::vector(row_size + 1, std::numeric_limits<int>::min() / 2);

    const auto [w_match, w_mismatch, w_open, w_extend] = params;
    for (auto i = 1; i < row_size; i++) {
      for (auto j = 1; j < col_size; j++) {
        // step diag
        const auto step_diag
          = score[i - 1][j - 1]
            + (ref[i - 1] == alt[j - 1] ? w_match : w_mismatch);

        // step down
        const auto gap_open_down = score[i - 1][j] + w_open;
        best_gap_down[j] += w_extend;
        if (gap_open_down > best_gap_down[j]) {
          best_gap_down[j] = gap_open_down;
          gap_size_down[j] = 1;
        } else
          gap_size_down[j]++;
        const auto step_down = best_gap_down[j];
        const auto step_down_size = gap_size_down[j];

        // step right
        const auto gap_open_right = score[i][j - 1] + w_open;
        best_gap_right[i] += w_extend;
        if (gap_open_right > best_gap_right[i]) {
          best_gap_right[i] = gap_open_right;
          gap_size_right[i] = 1;
        } else
          gap_size_right[i]++;
        const auto step_right = best_gap_right[i];
        const auto step_right_size = gap_size_right[i];

        // priority here will be step diagonal, step down, step right
        if (step_diag >= step_down && step_diag >= step_right) {
          score[i][j] = step_diag;
          trace[i][j] = 0;
        } else if (step_right >= step_down) {
          score[i][j] = step_right;
          trace[i][j] = -step_right_size;
        } else {
          score[i][j] = step_down;
          trace[i][j] = step_down_size;
        }
      }
    }
  }

  static auto
  calculate_cigar(std::vector<std::vector<int>>& score,
                  std::vector<std::vector<int>>& trace) {
    const auto ref_size = score.size() - 1;
    const auto alt_size = score.front().size() - 1;

    auto max_score = std::numeric_limits<int>::min();
    auto segment_len = 0;

    // look for the largest score on the rightmost column. we use >= combined
    // with the traversal direction to ensure that if two scores are equal, the
    // one closer to diagonal gets picked
    auto pos_i = 0;
    for (auto i = 1; i <= ref_size; i++) {
      const auto cur_score = score[i][alt_size];
      if (cur_score >= max_score) {
        max_score = cur_score;
        pos_i = i;
      }
    }

    // now look for a larger score on the bottom-most row
    auto pos_j = alt_size;
    auto diff = [](auto x, auto y) { return x > y ? x - y : y - x; };
    for (auto j = 1; j <= alt_size; j++) {
      const auto cur_score = score[ref_size][j];
      if (cur_score > max_score
          || (cur_score == max_score
              && diff(ref_size, j) < diff(pos_i, pos_j))) {
        max_score = cur_score;
        pos_i = ref_size;
        pos_j = j;
        // end of alternate is overhanging; we will just record it as 'M'
        // segment
        segment_len = alt_size - j;
      }
    }

    auto cigar = Cigar{};
    if (segment_len > 0) {
      cigar.emplace_back(segment_len, 'S');
      segment_len = 0;
    }

    auto state = 'M';
    do {
      const auto cur_trace = trace[pos_i][pos_j];
      const auto [new_state, step_size] = [cur_trace] {
        if (cur_trace > 0)
          return std::pair{'D', cur_trace};
        else if (cur_trace < 0)
          return std::pair{'I', -cur_trace};
        else
          return std::pair{'M', 1};
      }();

      // move to next best location in the sw matrix
      switch (new_state) {
        case 'M':
          pos_i--;
          pos_j--;
          break;
        case 'I':
          pos_j -= step_size;
          break;
        case 'D':
          pos_i -= step_size;
          break;
        default:
          break;
      }

      if (new_state == state)
        segment_len += step_size;
      else {
        cigar.emplace_back(segment_len, state);
        segment_len = step_size;
        state = new_state;
      }
    } while (pos_i > 0 && pos_j > 0);

    cigar.emplace_back(segment_len, state);
    auto align_offset = pos_i;
    if (pos_j > 0)
      cigar.emplace_back(pos_j, 'S');

    cigar.reverse();
    return std::pair{align_offset, cigar};
  }

 public:
  static auto
  align(std::string_view ref, std::string_view alt,
        Parameters params = NEW_SW_PARAMETERS) {
    assert(!ref.empty() && !alt.empty());

    if (alt.size() == ref.size() && well_match(ref, alt))
      return std::pair{0, Cigar(std::to_string(ref.size()) + 'M')};

    auto score = std::vector(ref.size() + 1, std::vector(alt.size() + 1, 0));
    auto trace = std::vector(ref.size() + 1, std::vector(alt.size() + 1, 0));
    calculate_matrix(ref, alt, score, trace, params);
    return calculate_cigar(score, trace);
  }
};

}  // namespace biovoltron
