#include <biovoltron/algo/sort/psais_sorter.hpp>
#include <biovoltron/utility/istring.hpp>
#include <experimental/random>
#include <catch.hpp>
#include <chrono>
#include <algorithm>

using namespace biovoltron;

TEST_CASE("PsaisSorter") {
  auto gen_dna_seq = [](int len) {
    auto seq = std::string{};
    while (len--)
      seq += "ACGT"[std::experimental::randint(0, 3)];
    return seq;
  };

  int len = std::experimental::randint(100'000, 200'000);
  const auto seq = gen_dna_seq(len);
  const auto seq_sv = std::string_view{seq};
  const auto ref = Codec::to_istring(seq);

  auto sa = PsaisSorter<>::get_sa(ref);
#pragma omp parallel for
  for (int i = 1; i < sa.size(); i++)
    REQUIRE(seq_sv.substr(sa[i - 1]) < seq_sv.substr(sa[i]));
}

TEST_CASE("PsaisSorter large testcase") {
  auto gen_dna_seq = [](int len) {
    auto seq = std::string{};
    while (len--)
      seq += "ACGT"[std::experimental::randint(0, 3)];
    return seq;
  };

  int len = std::experimental::randint(10'000'000, 20'000'000);
  const auto seq = gen_dna_seq(len);
  const auto seq_sv = std::string_view{seq};
  const auto ref = Codec::to_istring(seq);

  auto sa = PsaisSorter<>::get_sa(ref);
#pragma omp parallel for
  for (int i = 1; i < sa.size(); i++)
    REQUIRE(seq_sv.substr(sa[i - 1]) < seq_sv.substr(sa[i]));
}
