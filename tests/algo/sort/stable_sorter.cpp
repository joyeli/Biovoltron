#include <biovoltron/algo/sort/stable_sorter.hpp>
#include <biovoltron/utility/istring.hpp>
#include <catch.hpp>
#include <chrono>
#include <algorithm>
#include <iostream> 

using biovoltron::StableSorter;
using namespace biovoltron;
using namespace std::chrono;

TEST_CASE("StableSorter::get_sa - Sorts stably", "[StableSorter]") {
  SECTION("Sort all the way till the end of the sequence") {
    auto seq = Codec::to_istring("acgtaacca");
    auto sa = StableSorter<std::uint32_t>::get_sa(seq); 
    REQUIRE(std::ranges::is_sorted(sa, [seq](auto i, auto j) {
      return seq.substr(i) < seq.substr(j);
    }));
  }

  SECTION("Sort 2 base") {
    auto seq = Codec::to_istring("acgtaacca");
    auto sa = StableSorter<std::uint32_t>::get_sa(seq, 2);
    REQUIRE(std::ranges::is_sorted(sa, [seq](auto i, auto j) {
      return seq.substr(i, 2) < seq.substr(j, 2);
    }));
  }

  // SECTION("Performance test") {
  //   const auto length = 1024 * 1024; 
  //   auto seq = istring{}; seq.reserve(length);
  //   std::default_random_engine eng;
  //   std::uniform_int_distribution<std::int8_t> dist(0, 3); 
  //   std::generate_n(
  //     std::back_inserter(seq), 
  //     length, 
  //     [&eng, &dist]() { return dist(eng); });
  //
  //   auto start = high_resolution_clock::now();
  //   auto sa = Sorter::get_sa<std::uint32_t>(seq, 256);
  //   auto end = high_resolution_clock::now();
  //   auto dur = duration_cast<milliseconds>(end - start);
  //   std::cerr << "Get suffix array of sequence with length " << length 
  //             << ", elapsed time: " << dur.count() << " ms.\n";
  //   REQUIRE(std::ranges::is_sorted(sa, [seq](auto i, auto j) {
  //     return seq.substr(i, 256) < seq.substr(j, 256);
  //   }));
  // }
}
