#pragma once

#include <vector>
#include <array>
#include <ranges>
#include <concepts>
#include <memory>
#include <algorithm>
#include <utility>

#include <simdpp/simd.h>
#include <spdlog/spdlog.h>

#include <biovoltron/utility/istring.hpp>

namespace biovoltron::detail {

/* how many bytes does a simd vector have */
constexpr auto VECTOR_SIZE = static_cast<std::size_t>(SIMDPP_FAST_INT8_SIZE);
/* how many bases inside a simd vector */
constexpr auto BASE_IN_VECTOR = VECTOR_SIZE << 2; // 2 bit per base
/* hit count inside a byte */

constexpr auto match_in_byte = []() constexpr {
  constexpr auto sz = (1u << (sizeof(std::uint8_t) * 8));
  std::array<uint8_t, sz> hit { 0 };
  std::ranges::for_each(std::views::iota(0u, sz), [&hit](auto idx) constexpr {
    hit[idx] = hit[idx >> 2] + ((idx & 3) == 3);
  });
  return hit;
}();

using SIMD_VECTOR = simdpp::uint8<VECTOR_SIZE>;

/* mask some bases when comparing in specific length */
inline static std::vector<SIMD_VECTOR> erase_mask = []() {
  auto mask = std::vector<SIMD_VECTOR>(BASE_IN_VECTOR + 1);
  auto ptr = static_cast<uint8_t*>(
    std::aligned_alloc(VECTOR_SIZE, VECTOR_SIZE)
  );
  memset(ptr, 0xff, VECTOR_SIZE);
  mask[0] = simdpp::load(ptr);
  for (auto i = 1u; i <= BASE_IN_VECTOR; i++) {
    int idx = (i - 1) >> 2;
    *(ptr + idx) >>= 2;
    mask[i] = simdpp::load(ptr);
  }
  std::free(ptr);
  return mask;
}();

#ifdef DEBUG
template<std::size_t... IDX>
inline auto to_string(const SIMD_VECTOR& v, std::index_sequence<IDX...>) {
  std::string s;
  auto foo = [&s](uint8_t x) {
    for (int i = 6; i >= 0; i -= 2) {
      s += Codec::to_char((x >> i) & 3);
    }
  };
  ((foo(simdpp::extract<IDX>(v))), ...);
  return s;
}

inline auto to_string(const SIMD_VECTOR& v) {
  return to_string(v, std::make_index_sequence<VECTOR_SIZE>());
}
#endif

/**
 * @brief make a `SIMD_VECTOR` for given sequence
 *
 * @param seq The sequence
 * @param vec_size How many `SIMD_VECTOR` should be made. The value of -1
 * indicate no specify the size, will give the result of
 * `std::vector<SIMD_VECTOR>` which store all bases in `seq`
 * @return `std::vector<SIMD_VECTOR>`
 */
template<std::ranges::view SeqView>
auto make_simd_vector(SeqView seq, std::size_t vec_size = -1) {

  auto sz = seq.size();
  auto need_size = sz / BASE_IN_VECTOR + (sz % BASE_IN_VECTOR != 0);
  vec_size = std::min(vec_size, need_size);
  auto vecs = std::vector<SIMD_VECTOR> {};

  auto fill_byte = [](SeqView s, std::uint8_t* ptr) {
    *ptr = 0;
    for (auto i = 0u; i < 4; i++) {
      *ptr <<= 2;
      if (i < s.size()) {
        if constexpr (std::same_as<SeqView, istring_view>) {
          *ptr |= s[i];
        } else {
          *ptr |= Codec::to_int(s[i]);
        }
      }
    }
  };

  for (auto i = 0u; i < vec_size; i++) {
    auto ptr = static_cast<std::uint8_t*>(
      std::aligned_alloc(VECTOR_SIZE, VECTOR_SIZE)
    );
    for (auto k = 0u; k < VECTOR_SIZE; k++) {
      auto offset = i * BASE_IN_VECTOR + k * 4;
      if (offset >= sz) {
        break;
      }
      auto s = seq.substr(offset, 4);
      fill_byte(s, ptr + k);
    }
    auto v = SIMD_VECTOR(simdpp::load(ptr));
    vecs.push_back(v);
    std::free(ptr);
  }
  return vecs;
};

auto shift_left(const SIMD_VECTOR& v) -> SIMD_VECTOR {
  if constexpr (VECTOR_SIZE == 16) {
    return v << 2 | simdpp::move16_l<1>(v >> 6);
  } else {
    /* simdpp::move16_l will perform each 128 bits seperately */

    /* how many 128 bits should be move left */
    constexpr auto level = VECTOR_SIZE >> 4;

    /* store the 2 left most bits in v[8], v[16]... */
    auto r = SIMD_VECTOR { v << 2 | simdpp::move16_l<1>(v >> 6) };

    /* restore the 2 left most bits in the right most bist of v[7], v[15]...*/
    auto restore_bit = [&]<std::size_t... IDX>(std::index_sequence<IDX...>) {
      (
        [&]() {
        constexpr auto pos = static_cast<uint8_t>((IDX + 1) << 4);
        static_assert(0 < pos && pos < VECTOR_SIZE);
        auto a = simdpp::extract<pos - 1>(r);
        auto b = simdpp::extract<pos>(SIMD_VECTOR(v >> 6));
        r = simdpp::insert<pos - 1>(r, a | b);
        }(), ...  
      );
    };
    restore_bit(std::make_index_sequence<level - 1>());
    return r;
  }
}

/**
  * @brief calculate match bases between two `SIMD_VECTOR`
  *
  * @param v1 vector 1
  * @param v2 vector 2
  * @param cmp_len the base we used to calcuate the matched base
  * @return `std::size_t`
  */
inline auto cal_match(const SIMD_VECTOR& v1,
                      const SIMD_VECTOR& v2,
                      const std::size_t cmp_len) -> std::size_t {

  const auto& mask = erase_mask[cmp_len];
  auto v = SIMD_VECTOR(~((v1 | mask) ^ (v2 & (~mask))));
  auto counter = [&]<std::size_t... IDX>(std::index_sequence<IDX...>) {
    auto count = 0;
    ((count += match_in_byte[simdpp::extract<IDX>(v)]), ...);
    return count;
  };
  // ? benchmark: the simdpp builtin function and fold experssion
  // simdpp::for_each(v, [&](const uint8_t& u) {
  //   count += match_in_byte[u];
  // });
  return counter(std::make_index_sequence<VECTOR_SIZE>());
}

/**
 * @brief calculate the similarity between two `SIMD_VECTOR`
 *
 * @param v1 vector 1
 * @param v2 vector 2
 * @param cmp_len the base we used to calcuate the similarity
 * @return `double` the similarity
 */
inline auto cal_similarity(const SIMD_VECTOR& v1,
                           const SIMD_VECTOR& v2,
                           const std::size_t cmp_len) -> double {
  assert(0 < cmp_len && cmp_len <= BASE_IN_VECTOR);
  return static_cast<double>(cal_match(v1, v2, cmp_len)) / cmp_len;
}

}