#define SPDLOG_ACTIVE_LEVEL SPLLOG_LEVEL_DEBUG
#include <spdlog/sinks/ostream_sink.h>
#include <biovoltron/algo/sort/psais_sorter.hpp>
#include <biovoltron/algo/align/exact_match/fm_index.hpp>
#include <catch.hpp>
#include <experimental/random>
#include <iomanip>
#include <iostream>

using namespace biovoltron;

TEST_CASE("FMIndex") {
  spdlog::set_level(spdlog::level::info);
  
  auto gen_dna_seq = [](int len) {
    auto seq = std::string{};
    while (len--)
      seq += "ATGC"[std::experimental::randint(0, 3)];
    return seq;
  };

  const int LOOKUP_LEN = 8;
  uint32_t len = std::experimental::randint(500, 1000);
  auto seq = gen_dna_seq(len);

  // append LOOKUP_LEN 'A' after seq
  for (int i = 0; i < LOOKUP_LEN; i++)
    seq += 'A';

  const auto ref = Codec::to_istring(seq);
  auto fmidx = biovoltron::FMIndex{
    .LOOKUP_LEN = LOOKUP_LEN
  };
  fmidx.build(ref);

  SECTION("suffix array") {
    for (int i = 1; i < fmidx.sa_.size(); i++)
      REQUIRE(seq.substr(fmidx.sa_[i - 1]) < seq.substr(fmidx.sa_[i]));
  }

  SECTION("lookup table") {
    int num_suffixes = 0;
    for (int i = 0; i < (1 << (LOOKUP_LEN * 2)); i++) {
      auto sd = Codec::to_string(Codec::rhash(i, LOOKUP_LEN));
      const auto [beg, end, offset] = fmidx.get_range(Codec::to_istring(sd), -1);
      num_suffixes += end - beg;
      auto hits = fmidx.get_offsets(beg, end);
      for (auto &v : hits)
        REQUIRE(seq.substr(v, LOOKUP_LEN) == sd);
    }
    REQUIRE(num_suffixes + LOOKUP_LEN == fmidx.sa_.size());
  }

  SECTION("parallel occ table") {
    auto cnt1 = std::array<std::uint32_t, 4>{};
    auto cnt2 = std::array<std::uint8_t, 4>{};
    auto &[occ1, occ2] = fmidx.occ_;
    for (int i = 0; i < fmidx.sa_.size(); i++) {
      if (i % fmidx.OCC1_INTV == 0) {
        cnt2 = {};
        REQUIRE(occ1[i / fmidx.OCC1_INTV] == cnt1);
      }
      if (i % fmidx.OCC2_INTV == 0)
        REQUIRE(occ2[i / fmidx.OCC2_INTV] == cnt2);

      const auto sa_v = fmidx.sa_[i];
      if (sa_v != 0) {
        const auto c = ref[sa_v - 1];
        cnt1[c]++;
        cnt2[c]++;
      }
    }
  }

  SECTION("random query") {
    int q = 100;
    while (q--) {
      int seed_len = std::experimental::randint(5, 13);
      const auto seed_seq = gen_dna_seq(seed_len);
      const auto seed = Codec::to_istring(seed_seq);
      const auto [beg, end, offset] = fmidx.get_range(seed, 0);
      auto hits = fmidx.get_offsets(beg, end);

      // check hits with brute-force
      int num_hits = 0;
      for (int i = 0; i < (int)seq.size(); i++)
        if (seq.substr(i, seed_len) == seed_seq)
          num_hits++;

      REQUIRE(hits.size() == num_hits);
      for (const auto &hit : hits)
        REQUIRE(seq.substr(hit, seed_len) == seed_seq);
    }
  }
}

TEST_CASE("FMIndex with SA sampling") {
  auto gen_dna_seq = [](int len) -> std::string {
    auto seq = std::string{};
    while (len--)
      seq += "ATGC"[std::experimental::randint(0, 3)];
    return seq;
  };

  const int LOOKUP_LEN = 8;
  uint32_t len = std::experimental::randint(500, 1000);
  auto seq = gen_dna_seq(len);

  // append LOOKUP_LEN 'A' after seq
  for (int i = 0; i < LOOKUP_LEN; i++)
    seq += 'A';

  const auto ref = Codec::to_istring(seq);
  auto fmidx = biovoltron::FMIndex<8>{
    .LOOKUP_LEN = LOOKUP_LEN
  };

  fmidx.build(ref);

  int q = 100;
  while (q--) {
    int seed_len = std::experimental::randint(5, 18);
    auto seed_seq = gen_dna_seq(seed_len);

    const auto seed = Codec::to_istring(seed_seq);
    auto [beg, end, offs] = fmidx.get_range(seed, 0);
    auto hits = fmidx.get_offsets(beg, end);

    // check hits with brute-force
    uint32_t checksum = 0, num_hits = 0;
#pragma omp parallel for reduction(+:checksum, num_hits)
    for (uint32_t i = 0; i < (uint32_t)seq.size(); i++) {
      if (seq.substr(i, seed_len) == seed_seq) {
        checksum += i;
        num_hits++;
      }
    }

    REQUIRE(hits.size() == num_hits);
    uint32_t sum = 0;
    for (const auto &hit : hits) {
      REQUIRE(seq.substr(hit, seed_len) == seed_seq);
      sum += hit;
    }
    REQUIRE(checksum == sum);
  }
}

TEST_CASE("FMIndex value sample") {
  auto gen_dna_seq = [](int len) -> std::string {
    auto seq = std::string{};
    while (len--)
      seq += "ATGC"[std::experimental::randint(0, 3)];
    return seq;
  };

  constexpr int SA_INTV = 8;
  const int LOOKUP_LEN = 0;
  uint32_t len = std::experimental::randint(500, 1000);
  auto seq = gen_dna_seq(len);

  // append LOOKUP_LEN 'A' after seq
  for (int i = 0; i < LOOKUP_LEN * 2; i++)
    seq += 'A';

  const auto ref = Codec::to_istring(seq);
  auto fmidx = biovoltron::FMIndex<SA_INTV>{
    .LOOKUP_LEN = LOOKUP_LEN
  };

  auto ori_sa = PsaisSorter<uint32_t>::get_sa(ref, istring::npos);
  fmidx.build(ref, ori_sa);

  SECTION("suffix array") {
    for (int i = 1; i < fmidx.sa_.size(); i++)
      REQUIRE(seq.substr(fmidx.sa_[i - 1]) < seq.substr(fmidx.sa_[i]));
  }

  SECTION("lookup table") {
    int num_suffixes = 0;
    for (int i = 0; i < (1 << (LOOKUP_LEN * 2)); i++) {
      auto sd = Codec::to_string(Codec::rhash(i, LOOKUP_LEN));
      const auto [beg, end, offset] = fmidx.get_range(Codec::to_istring(sd), 0);
      num_suffixes += end - beg;
      auto hits = fmidx.get_offsets(beg, end);
      for (auto &v : hits)
        REQUIRE(seq.substr(v, LOOKUP_LEN) == sd);
    }
    REQUIRE(num_suffixes + LOOKUP_LEN == fmidx.bwt_.size());
  }

  SECTION("occ table") {
    auto cnt1 = std::array<std::uint32_t, 4>{};
    auto cnt2 = std::array<std::uint8_t, 4>{};
    auto &[occ1, occ2] = fmidx.occ_;
    for (int i = 0; i < ori_sa.size(); i++) {
      if (i % fmidx.OCC1_INTV == 0) {
        cnt2 = {};
        REQUIRE(occ1[i / fmidx.OCC1_INTV] == cnt1);
      }
      if (i % fmidx.OCC2_INTV == 0)
        REQUIRE(occ2[i / fmidx.OCC2_INTV] == cnt2);

      const auto sa_v = ori_sa[i];
      if (sa_v != 0) {
        const auto c = ref[sa_v - 1];
        cnt1[c]++;
        cnt2[c]++;
      }
    }
  }

  SECTION("b") {
    REQUIRE(ori_sa.size() == fmidx.b_.size());
    for (int i = 0; i < fmidx.b_.size(); i++) {
      REQUIRE(fmidx.b_[i] == (ori_sa[i] % SA_INTV == 0));
    }
    auto cnt = uint32_t{};
    for (int i = 0; i < fmidx.b_.size(); i++) {
      if (i % fmidx.B_OCC_INTV == 0)
        REQUIRE(cnt == fmidx.b_occ_[i / fmidx.B_OCC_INTV]);
      cnt += fmidx.b_[i];
    }
  }
}
