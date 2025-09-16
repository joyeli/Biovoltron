#include <biovoltron/algo/sort/kiss_sorter/kiss2_sorter.hpp>
#include <biovoltron/utility/istring.hpp>
#include <experimental/random>
#include <catch.hpp>
// #include <catch_amalgamated.hpp>
#include <algorithm>
#include <vector>
#include <mutex>

using namespace biovoltron;

TEST_CASE("Kiss2Sorter DNA") {
  auto gen_dna_seq = [](int len) {
    auto seq = std::string{};
    while (len--)
      seq += "ACGT"[std::experimental::randint(0, 3)];
    return seq;
  };

  int len = std::experimental::randint(100'000, 200'000);
  auto k = size_t{256};
  const auto seq = gen_dna_seq(len);
  const auto seq_sv = std::string_view{seq};
  const auto ref = Codec::to_istring(seq);

  auto sa = KISS2Sorter<>::get_sa(ref, k, 24);
  std::vector<int> failed_indices;
  std::mutex vector_mutex;
#pragma omp parallel for
  for (int i = 1; i < sa.size(); i++) {
    if (!(seq_sv.substr(sa[i - 1], k) <= seq_sv.substr(sa[i], k))) {
        std::lock_guard<std::mutex> lock(vector_mutex);
        failed_indices.push_back(i);
    }
  }
  if (!failed_indices.empty()) {
    INFO("Test failed at indices: " << Catch::Detail::stringify(failed_indices));
  }
  REQUIRE(failed_indices.empty());
}

TEST_CASE("Kiss2Sorter DNA large") {
  auto gen_dna_seq = [](int len) {
    auto seq = std::string{};
    while (len--)
      seq += "ACGT"[std::experimental::randint(0, 3)];
    return seq;
  };

  int len = std::experimental::randint(10'000'000, 20'000'000);
  auto k = size_t{256};
  const auto seq = gen_dna_seq(len);
  const auto seq_sv = std::string_view{seq};
  const auto ref = Codec::to_istring(seq);

  auto sa = KISS2Sorter<>::get_sa(ref, k, 24);
  std::vector<int> failed_indices;
  std::mutex vector_mutex;
#pragma omp parallel for
  for (int i = 1; i < sa.size(); i++) {
    if (!(seq_sv.substr(sa[i - 1], k) <= seq_sv.substr(sa[i], k))) {
        std::lock_guard<std::mutex> lock(vector_mutex);
        failed_indices.push_back(i);
    }
  }
  if (!failed_indices.empty()) {
    INFO("Test failed at indices: " << Catch::Detail::stringify(failed_indices));
  }
  REQUIRE(failed_indices.empty());
}

TEST_CASE("Kiss2Sorter general") {
  auto gen_seq = [](int len) {
    auto seq = std::string{};
    while (len--)
      seq += std::experimental::randint(65, 68); // 'A' ~ 'D'
    return seq;
  };

  int len = std::experimental::randint(100'000, 200'000);
  auto k = size_t{256};
  const auto seq = gen_seq(len);
  const auto seq_sv = std::string_view{seq};

  auto sa = KISS2Sorter<>::get_suffix_array(seq_sv, k);
  std::vector<int> failed_indices;
  std::mutex vector_mutex;
#pragma omp parallel for
  for (int i = 1; i < sa.size(); i++) {
    if (!(seq_sv.substr(sa[i - 1], k) <= seq_sv.substr(sa[i], k))) {
        std::lock_guard<std::mutex> lock(vector_mutex);
        failed_indices.push_back(i);
    }
  }
  if (!failed_indices.empty()) {
    INFO("Test failed at indices: " << Catch::Detail::stringify(failed_indices));
  }
  REQUIRE(failed_indices.empty());
}

TEST_CASE("Kiss2Sorter general large") {
  auto gen_seq = [](int len) {
    auto seq = std::string{};
    while (len--)
      seq += std::experimental::randint(65, 68); // 'A' ~ 'D'
    return seq;
  };

  int len = std::experimental::randint(10'000'000, 20'000'000);
  auto k = size_t{256};
  const auto seq = gen_seq(len);
  const auto seq_sv = std::string_view{seq};

  auto sa = KISS2Sorter<>::get_suffix_array(seq_sv, k);
  std::vector<int> failed_indices;
  std::mutex vector_mutex;
#pragma omp parallel for
  for (int i = 1; i < sa.size(); i++) {
    if (!(seq_sv.substr(sa[i - 1], k) <= seq_sv.substr(sa[i], k))) {
        std::lock_guard<std::mutex> lock(vector_mutex);
        failed_indices.push_back(i);
    }
  }
  if (!failed_indices.empty()) {
    INFO("Test failed at indices: " << Catch::Detail::stringify(failed_indices));
  }
  REQUIRE(failed_indices.empty());
}