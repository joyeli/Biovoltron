#pragma once

#include <intel_native/smithwaterman/avx2-smithwaterman.h>

#include <biovoltron/algo/align/inexact_match/smithwaterman.hpp>
#include <memory>

namespace biovoltron {

/**
 * @ingroup align
 *
 * @brief TBA
 */
struct AvxSmithWaterman {
  static auto
  align(std::string_view ref, std::string_view alt,
        SmithWaterman::Parameters params = SmithWaterman::NEW_SW_PARAMETERS) {
    assert(!ref.empty() && !alt.empty());

    if (alt.size() == ref.size() && SmithWaterman::well_match(ref, alt))
      return std::pair{0, Cigar(std::to_string(ref.size()) + 'M')};

    const auto [match, mismatch, open, extend] = params;
    const auto refLength = ref.length();
    const auto altLength = alt.length();
    auto cigarArray
      = std::make_unique<char[]>(2 * std::max(refLength, altLength));
    auto count = 0;
    return std::pair{
      runSWOnePairBT_avx2(match, mismatch, open, extend, (uint8_t*)ref.data(),
                          (uint8_t*)alt.data(), refLength, altLength, 9,
                          cigarArray.get(), (int16_t*)&count),
      Cigar(cigarArray.get())};
  }
};

}  // namespace biovoltron
