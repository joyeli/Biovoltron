#pragma once

#include <ssw.h>

#include <array>
#include <string>

namespace biovoltron {

/**
 * @ingroup align
 *
 * @brief SSE-accelerated Smith–Waterman wrapper using the SSW library.
 *
 * @details
 * This struct provides a thin C++ wrapper over the SSW (Striped Smith–Waterman) library
 * to perform local alignment with affine gap penalties. It exposes:
 *  - a static scoring matrix (match/mismatch/ambiguous),
 *  - helpers to create an SSW profile from a read sequence,
 *  - an alignment routine that returns key coordinates and a CIGAR string.
 *
 * The expected alphabet is A/C/G/T (encoded as 0..3) plus an ambiguous base (index 4).
 * The scoring matrix is 5x5 and constructed according to the constants below.
 */
struct SseSmithWaterman {
  /// Match reward.
  constexpr static auto w_match = 1;
  /// Mismatch penalty (positive here; used as -w_mismatch in the matrix).
  constexpr static auto w_mismatch = 4;
  /// Gap-open penalty.
  constexpr static auto w_open = 6;
  /// Gap-extension penalty.
  constexpr static auto w_extend = 1;
  /// Ambiguous base penalty (applied whenever a fifth symbol is involved).
  constexpr static auto w_ambig = 1;

  /**
   * @brief 5x5 substitution matrix for A/C/G/T/N (or ambiguous) alphabet.
   *
   * @details
   * Layout is row-major. For bases 0..3 (A,C,G,T) the diagonal gets +w_match,
   * off-diagonals get -w_mismatch. Any cell involving the 5th symbol (index 4)
   * receives -w_ambig. The SSW library expects an int8 matrix.
   */
  constexpr static auto mat = [] {
    auto mat = std::array<std::int8_t, 25>{};
    auto k = 0;
    for (auto i = 0; i < 4; i++) {
      for (auto j = 0; j < 4; j++) mat[k++] = i == j ? w_match : -w_mismatch;
      mat[k++] = -w_ambig;
    }
    for (auto i = 0; i < 5; i++) mat[k++] = -w_ambig;
    return mat;
  }();

  /**
   * @brief Result of a single SW alignment.
   *
   * @param score  Best local alignment score.
   * @param score2 Secondary best score (as returned by SSW).
   * @param ref_beg 0-based start on the reference of the best alignment.
   * @param ref_end 0-based end on the reference (inclusive).
   * @param read_beg 0-based start on the read of the best alignment.
   * @param read_end 0-based end on the read (inclusive).
   * @param ref_end2 Reference end position for the suboptimal alignment.
   * @param cigar CIGAR string including soft clips for overhangs (if any).
   */
  struct SWResult {
    int score{};
    int score2{};
    int ref_beg{};
    int ref_end{};
    int read_beg{};
    int read_end{};
    int ref_end2{};
    std::string cigar;
  };

 private:
  /**
   * @brief Convert an SSW raw result into SWResult and destroy the SSW object.
   *
   * @param res       Pointer returned by `ssw_align`.
   * @param read_size Length of the query/read (to compute end soft-clip).
   * @return Parsed SWResult with positions and CIGAR string.
   *
   * @details
   * If `read_beg1 > 0`, a leading soft clip is emitted. The body CIGAR is
   * reconstructed from SSW's internal (len,op) encoding, and a trailing soft
   * clip is emitted if `read_end1` does not reach `read_size - 1`.
   * The function takes ownership of `res` and calls `align_destroy(res)`.
   */
  static auto
  extract_result(s_align* res, int read_size) {
    const auto [score1, score2, ref_beg1, ref_end1, read_beg1, read_end1,
                ref_end2, icigar, cigar_size]
      = *res;
    auto result = SWResult{score1,    score2,    ref_beg1, ref_end1,
                           read_beg1, read_end1, ref_end2, {}};
    if (auto& cigar = result.cigar; cigar_size > 0) {
      if (read_beg1 > 0)
        cigar += std::to_string(read_beg1) + 'S';

      for (auto i = 0; i < cigar_size; i++)
        cigar += std::to_string(cigar_int_to_len(icigar[i]))
                 + cigar_int_to_op(icigar[i]);

      if (const auto end = read_size - read_end1 - 1; end > 0)
        cigar += std::to_string(end) + 'S';
    }
    align_destroy(res);
    return result;
  }

 public:
  /**
   * @brief Build an SSW profile for a given read/query.
   *
   * @tparam R A contiguous range whose value_type is `std::int8_t`
   *           (e.g., `std::vector<std::int8_t>` or `std::span<const std::int8_t>`).
   * @param read Encoded read sequence (A/C/G/T as 0..3, ambiguous as 4).
   * @return `s_profile*` created by `ssw_init`; caller passes it to `align`.
   *
   * @note The caller does not own the matrix; it uses the static `mat` defined
   *       in this struct. The returned profile must remain valid during `align`.
   */
  template<std::ranges::contiguous_range R>
    requires std::same_as<std::ranges::range_value_t<R>, std::int8_t>
  static auto
  get_profile(const R& read) {
    return ssw_init(read.data(), read.size(), mat.data(), 5, 0);
  }

  /**
   * @brief Run SSE-accelerated local alignment against a reference.
   *
   * @tparam R A contiguous range whose value_type is `std::int8_t`.
   * @param profile      SSW profile created by `get_profile`.
   * @param ref          Encoded reference sequence.
   * @param report_beg   If true, request begin positions from SSW. (default = true)
   * @param report_cigar If true, request the CIGAR string from SSW. (default = true)
   * @param min_score    Minimum score threshold to report an alignment. (default = 0)
   * @return `SWResult` with alignment score, coordinates, and CIGAR string.
   *
   * @details
   * Internally sets the SSW `flag` bits according to `report_beg`/`report_cigar`,
   * invokes `ssw_align`, then converts the raw result via `extract_result`.
   * Gap penalties are taken from the constants `w_open` and `w_extend`.
   */
  template<std::ranges::contiguous_range R>
    requires std::same_as<std::ranges::range_value_t<R>, std::int8_t>
  static auto
  align(const s_profile& profile, const R& ref, bool report_beg = true,
        bool report_cigar = true, int min_score = 0) {
    auto flag = 0;
    if (report_beg)
      flag |= 0x08;
    if (report_cigar)
      flag |= 0x0F;
    auto res = ssw_align(&profile, ref.data(), ref.size(), w_open, w_extend,
                         flag, min_score, 32767, profile.readLen / 2);
    return extract_result(res, profile.readLen);
  }
};

}  // namespace biovoltron
