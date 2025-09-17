#include <biovoltron/algo/align/inexact_match/smithwaterman_sse.hpp>
#include <biovoltron/utility/istring.hpp>
#include <catch.hpp>
#include <chrono>
#include <iostream>

using namespace biovoltron;

// NOTE: Assume alignment score < 255 (ssw_init: score_size = 0).
//  - S: soft clipping
//  - M: match / mismatch
//  - I: insertion
//  - D: deletion
TEST_CASE("SseSmithWaterman::align - Performs Smith-Waterman alignment using SSE", "[SseSmithWaterman]")
{
  auto ref = Codec::to_istring(
    "AATCGAAGGTCGTAAGGACACGGTTGAGCGTTCAGCGTTCATGTGAGTCCTCACCACTTATGGCTCCATAGCCT"
    "GCTATTTAAGTGGGTTACCGGTCTCCGCCAAGTAGCTGGTGTAAGAACACAGTAACTGAGCCCAGTGTGATCAG"
    "CCCTAACGAGGTAC");

  auto read = ref;

  SECTION("Same reads")
  {
    //auto start = std::chrono::high_resolution_clock::now();

    auto profile = SseSmithWaterman::get_profile(read);
    const auto result = SseSmithWaterman::align(profile, ref);

    CHECK(result.cigar == "162M");

    //auto end = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    //std::cout << "SmithWaterman_sse_same time: " << duration << " us" << std::endl;
  }

  SECTION("Substitutions")
  {
    //auto start = std::chrono::high_resolution_clock::now();

    for (auto idx = 70; idx < 80; ++idx) {
      ref[idx] = Codec::to_int('A');
      read[idx] = Codec::to_int('T');
    }

    auto profile = SseSmithWaterman::get_profile(read);
    const auto result = SseSmithWaterman::align(profile, ref);

    CHECK(result.cigar == "70M10I10D82M");

    //auto end = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    //std::cout << "SmithWaterman_sse_substi time: " << duration << " us" << std::endl;
  }

  SECTION("Deletion")
  {
    //auto start = std::chrono::high_resolution_clock::now();

    read.erase(read.begin() + 70);

    auto profile = SseSmithWaterman::get_profile(read);
    const auto result = SseSmithWaterman::align(profile, ref);

    CHECK(result.cigar == "70M1D91M");

    //auto end = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    //std::cout << "SmithWaterman_sse_del time: " << duration << " us" << std::endl;
  }

  SECTION("Insertion")
  {
    //auto start = std::chrono::high_resolution_clock::now();

    read.insert(read.begin() + 70, Codec::to_int('T'));

    auto profile = SseSmithWaterman::get_profile(read);
    const auto result = SseSmithWaterman::align(profile, ref);

    CHECK(result.cigar == "70M1I92M");

    //auto end = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    //std::cout << "SmithWaterman_sse_ins time: " << duration << " us" << std::endl;
  }

  SECTION("Mix")
  {
    //auto start = std::chrono::high_resolution_clock::now();
    
    for (auto idx = 0; idx < static_cast<int>(ref.size()); ++idx) {
      if ((idx > 10 && idx < 20) || (idx > 70 && idx < 80) ||
          (idx > 120 && idx < 130)) {
        ref[idx] = Codec::to_int('A');
        read[idx] = Codec::to_int('T');
      }
    }

    read.erase(read.begin() + 60);
    read.insert(read.begin() + 90, Codec::to_int('T'));

    auto profile = SseSmithWaterman::get_profile(read);
    const auto result = SseSmithWaterman::align(profile, ref);

    CHECK(result.cigar == "20S40M1D10M9I9D11M1I30M9I9D32M");

    //auto end = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    //std::cout << "SmithWaterman_sse_mix time: " << duration << " us" << std::endl;
  }
}
