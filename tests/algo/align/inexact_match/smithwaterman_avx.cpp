#include <biovoltron/algo/align/inexact_match/smithwaterman_avx.hpp>
#include <biovoltron/utility/istring.hpp>
#include <catch.hpp>
#include <chrono>
#include <iostream>

using namespace biovoltron;

TEST_CASE("Smith Waterman AVX")
{
  // Reference DNA sequence (length = 162)
  std::string ref =
    "AATCGAAGGTCGTAAGGACACGGTTGAGCGTTCAGCGTTCATGTGAGTCCTCACCACTTATGGCTCCATAGCCT"
    "GCTATTTAAGTGGGTTACCGGTCTCCGCCAAGTAGCTGGTGTAAGAACACAGTAACTGAGCCCAGTGTGATCAG"
    "CCCTAACGAGGTAC";

  // alt is a copy of ref used as query/read sequence
  auto alt = ref;

  SECTION("Same reads")
  {
    //auto start = std::chrono::high_resolution_clock::now();  // Start timer
    // Case: Identical sequences
    const auto [offset, cigar] = AvxSmithWaterman::align(ref, alt);

    REQUIRE(offset == 0);                    // Expect alignment from start
    REQUIRE(std::string(cigar) == "162M");   // Perfect match of 162 bases

    // Stop timer and print elapsed time
    //auto end = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
  
    //std::cout << "SmithWaterman_avx＿same time: " << duration << " us" << std::endl;
  }

  SECTION("Substitutions")
  {
    //auto start = std::chrono::high_resolution_clock::now();  // Start timer
    // Introduce mismatches between positions 70–79
    for (auto idx = 70; idx < 80; ++idx) {
      ref[idx] = 'A';
      alt[idx] = 'T';
    }

    const auto [offset, cigar] = AvxSmithWaterman::align(ref, alt);

    REQUIRE(offset == 0);
    REQUIRE(std::string(cigar) == "69M10D1M10I82M"); // Mismatches mapped to 10D + 10I

    // Stop timer and print elapsed time
    //auto end = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
  
    //std::cout << "SmithWaterman_avx_substi time: " << duration << " us" << std::endl;
  }

  SECTION("Deletion")
  {
    //auto start = std::chrono::high_resolution_clock::now();  // Start timer
    // Simulate deletion at position 70 in the alt sequence
    alt.erase(alt.begin() + 70);

    const auto [offset, cigar] = AvxSmithWaterman::align(ref, alt);

    REQUIRE(offset == 0);
    REQUIRE(std::string(cigar) == "70M1D91M");

    // Stop timer and print elapsed time
    //auto end = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
  
    //std::cout << "SmithWaterman_avx_del time: " << duration << " us" << std::endl;
  }

  SECTION("Insertion")
  {
    //auto start = std::chrono::high_resolution_clock::now();  // Start timer
    // Simulate insertion of 'T' at position 70 in the alt sequence
    alt.insert(alt.begin() + 70, 'T');

    const auto [offset, cigar] = AvxSmithWaterman::align(ref, alt);  // 🔧 FIX: use AVX version

    REQUIRE(offset == 0);
    REQUIRE(std::string(cigar) == "70M1I92M");

    // Stop timer and print elapsed time
    //auto end = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
  
    //std::cout << "SmithWaterman_avx_ins time: " << duration << " us" << std::endl;
  }

  SECTION("Mix")
  {
    //auto start = std::chrono::high_resolution_clock::now();  // Start timer
    // Introduce substitutions at 3 regions
    for (auto idx = 0; idx < static_cast<int>(ref.size()); ++idx) {
      if ((idx > 10 && idx < 20) || (idx > 70 && idx < 80) ||
          (idx > 120 && idx < 130)) {
        ref[idx] = 'A';
        alt[idx] = 'T';
      }
    }

    // Simulate deletion and insertion
    alt.erase(alt.begin() + 60);        // Deletion
    alt.insert(alt.begin() + 90, 'T');  // Insertion

    const auto [offset, cigar] = AvxSmithWaterman::align(ref, alt);  //FIX: use AVX version

    REQUIRE(offset == 0);
    REQUIRE(std::string(cigar) == "11M9D9I40M1D10M9D9I11M1I28M9D2M9I32M");
    // Complex alignment: soft clipping, matches, indels

    // Stop timer and print elapsed time
    //auto end = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
  
    //std::cout << "SmithWaterman_avx_mix time: " << duration << " us" << std::endl;
  }
}
