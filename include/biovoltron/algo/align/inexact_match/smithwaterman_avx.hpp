#pragma once

#include <intel_native/smithwaterman/avx2-smithwaterman.h>

#include <biovoltron/algo/align/inexact_match/smithwaterman.hpp>
#include <memory>

namespace biovoltron {

/**
 * @ingroup align
 *
 * @brief AVX2-accelerated Smith–Waterman wrapper.
 *
 * @details
 * This struct provides a C++ wrapper for an AVX2-optimized Smith–Waterman alignment
 * implementation (`runSWOnePairBT_avx2`). It performs local alignment with affine gap penalties
 * between two sequences and returns both the alignment offset and CIGAR string.
 *
 * It uses the same scoring parameters (`SmithWaterman::Parameters`) as the scalar
 * SmithWaterman implementation in this project but offloads computation to AVX2 SIMD intrinsics
 * for faster execution.
 */
struct AvxSmithWaterman {
  /**
   * @brief Align two sequences using AVX2-accelerated Smith–Waterman.
   *
   * @param ref    Reference sequence (as a `std::string_view`).
   * @param alt    Alternate/query sequence (as a `std::string_view`).
   * @param params Scoring parameters (match, mismatch, gap open, gap extend).
   *               Defaults to `SmithWaterman::NEW_SW_PARAMETERS`.
   *
   * @return A `std::pair<int, Cigar>` where:
   *         - First element: alignment offset (position in the reference where alignment starts).
   *         - Second element: CIGAR string describing the alignment.
   *
   * @details
   * - If the reference and alternate are the same length and differ by no more than
   *   `SmithWaterman::MAX_MISMATCHES`, the result is returned directly as a perfect match (`M`).
   * - Otherwise, the AVX2 routine `runSWOnePairBT_avx2` is called with:
   *   - Match score: `params.w_match`
   *   - Mismatch penalty: `params.w_mismatch`
   *   - Gap open penalty: `params.w_open`
   *   - Gap extend penalty: `params.w_extend`
   *   - Input sequences as raw `uint8_t*`
   *   - `cigarArray` buffer large enough for up to `2 * max(refLength, altLength)` operations
   * - The resulting CIGAR is built from the `cigarArray` buffer.
   *
   * @note
   * Requires CPU support for AVX2 instructions.
   */
  static auto
  align(std::string_view ref, std::string_view alt,
        SmithWaterman::Parameters params = SmithWaterman::NEW_SW_PARAMETERS) {
    assert(!ref.empty() && !alt.empty());

    // Fast path for nearly perfect matches
    if (alt.size() == ref.size() && SmithWaterman::well_match(ref, alt))
      return std::pair{0, Cigar(std::to_string(ref.size()) + 'M')};

    const auto [match, mismatch, open, extend] = params;
    const auto refLength = ref.length();
    const auto altLength = alt.length();

    // Allocate buffer for CIGAR output
    auto cigarArray
      = std::make_unique<char[]>(2 * std::max(refLength, altLength));
    auto count = 0;

    return std::pair{
      runSWOnePairBT_avx2(match, mismatch, open, extend,
                          (uint8_t*)ref.data(), (uint8_t*)alt.data(),
                          refLength, altLength, 9,
                          cigarArray.get(), (int16_t*)&count),
      Cigar(cigarArray.get())};
  }
};

}  // namespace biovoltron
