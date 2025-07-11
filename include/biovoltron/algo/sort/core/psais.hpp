#pragma once

#include "omp.h"

#include <biovoltron/container/xbit_vector.hpp>
#include <biovoltron/utility/thread_pool.hpp>
#include <biovoltron/utility/istring.hpp>

#include <ranges>
#include <numeric>
#include <tuple>
#include <vector>
#include <thread>
#include <future>

namespace biovoltron {

namespace psais {
  
#define NUM_THREADS        std::thread::hardware_concurrency()
#define INDUCE_NUM_THREADS std::min(NUM_THREADS, 16u)

  template<typename T>
  constexpr auto EMPTY = std::numeric_limits<T>::max();
  constexpr auto BLOCK_SIZE = 1u << 20;

  enum SUFFIX_TYPE { L_TYPE = 0, S_TYPE = 1 };
  using TypeVector
    = biovoltron::detail::XbitVector<1, std::uint8_t,
                                     std::allocator<std::uint8_t>>;

  template <typename size_type>
  auto
  is_LMS(const auto& T, size_type i) {
    auto n = (size_type)T.size();
    return i == n
           or (i > 0 and T[i - 1] == SUFFIX_TYPE::L_TYPE
               and T[i] == SUFFIX_TYPE::S_TYPE);
  }

  template <typename size_type>
  auto
  num_lms(const auto &T) {
    auto n = (size_type)T.size();
    auto len = size_type{};
#pragma omp parallel for num_threads(NUM_THREADS) reduction(+:len)
    for (auto i = size_type{}; i < n; i++)
      if (is_LMS(T, i))
        len++;
    return len;
  }

  auto
  get_type_block_range(auto num_items, auto num_blocks, auto block_idx) {
    // in order to make TypeVector(a.k.a XBitVector) thread-safe,
    // the element should be group by 8 bit at least
    auto N = (num_items + 7) / 8, m = num_blocks;
    auto counts = N / m, remain = N % m;
    auto L = counts * block_idx + std::min(block_idx, remain);
    auto R = L + counts + (block_idx < remain);
    return std::tuple{std::min(num_items, L * 8), std::min(num_items, R * 8)};
  }

  template <typename size_type>
  auto
  get_type_per_block(const std::ranges::random_access_range auto& S, auto& T) {
    auto n = (size_type)S.size();
    auto suffix_len = std::vector<size_type>(NUM_THREADS, 0);
#pragma omp parallel for num_threads(NUM_THREADS)
    for (auto tid = size_type{}; tid < NUM_THREADS; tid++) {
      auto [L, R] = get_type_block_range(n, NUM_THREADS, tid);
      if (L == R)
        continue;

      for (auto i = R - L - 1 - (R == n), same = 1; ~i; i--) {
        auto x = L + i;
        auto c1 = S[x], c2 = S[x + 1];
        T[x] = c1 == c2 ? T[x + 1] :
               c1 < c2  ? SUFFIX_TYPE::S_TYPE :
                          SUFFIX_TYPE::L_TYPE;

        if (same) {
          same &= (S[x] == S[x + 1]);
          suffix_len[tid] += same;
        }
      }
    }
    return suffix_len;
  }

  template <typename size_type>
  auto
  get_type_check_flip(const std::ranges::random_access_range auto& S, auto& T,
                      const auto& suffix_len) {
    auto n = (size_type)S.size();
    auto flip = std::vector<uint8_t>(NUM_THREADS, false);
    for (auto tid = size_type{NUM_THREADS - 2}; ~tid; tid--) {
      auto [L, R] = get_type_block_range(n, NUM_THREADS, tid + 1);
      if (L == R)
        continue;

      //        tid   tid + 1
      // ...-|-------|-------|-...
      //           x1 x2
      auto x1 = L - 1, x2 = L;
      if (S[x1] != S[x2])
        continue;

      auto should_flip = suffix_len[tid + 1] == R - L and flip[tid + 1];
      flip[tid] = T[x1] != (T[x2] ^ should_flip);
    }
    return flip;
  }

  template <typename size_type>
  void
  get_type_flip_block(const std::ranges::random_access_range auto& S, auto& T,
                      const auto& suffix_len, const auto& flip) {
    auto n = (size_type)S.size();
#pragma omp parallel for num_threads(NUM_THREADS)
    for (auto tid = size_type{}; tid < NUM_THREADS; tid++) {
      if (flip[tid]) {
        auto [L, R] = get_type_block_range(n, NUM_THREADS, tid);
        for (auto i = size_type{}; i < suffix_len[tid]; i++)
          T[R - i - 1] = not T[R - i - 1];
      }
    }
  }

  template <typename size_type>
  void
  get_type(const std::ranges::random_access_range auto& S, auto& T) {
    auto n = (size_type)S.size();

    // if the string S is empty, return an empty TypeVector
    if (n == 0)
      return;

    // trivial case, the end character is L_TYPE
    T[n - 1] = SUFFIX_TYPE::L_TYPE;

    // calculate suffix type per block independently
    // assume that the last character in the block is L_TYPE
    auto suffix_len = get_type_per_block<size_type>(S, T);

    // check the last character is S_TYPE or L_TYPE
    auto flip = get_type_check_flip<size_type>(S, T, suffix_len);

    // if the last character is S_TYPE in fact,
    // flip suffix in block with same character
    get_type_flip_block<size_type>(S, T, suffix_len, flip);
  }

  template<auto induce_type, typename size_type>
  auto
  get_bucket(const std::vector<size_type>& BA_) {
    if (BA_.size() == 0)
      return BA_;

    if constexpr (induce_type == SUFFIX_TYPE::L_TYPE) {
      auto BA = std::vector<size_type>(BA_.size(), 1);
#pragma omp parallel for num_threads(NUM_THREADS)
      for (auto i = size_type{1}; i < BA_.size(); i++) BA[i] = BA_[i - 1];
      return BA;
    } else {
      return BA_;
    }
  }

  auto
  split_into_chunks(auto num_items, auto mem_size,
                   auto num_items_per_chunk) {
    
    auto num_chunks = (mem_size - 1) / num_items_per_chunk + 1;
    auto chunk_size = (num_items - 1) / num_chunks + 1;
    auto num_threads = num_chunks < NUM_THREADS ? num_chunks : NUM_THREADS;
    return std::tuple{num_chunks, chunk_size, num_threads};
  }

  template<typename size_type, typename F>
  auto
  get_local_bucket(const std::ranges::random_access_range auto& S, size_type K,
                   F&& check) {
    auto n = (size_type)S.size();
    auto [num_chunks, chunk_size, num_threads]
      = split_into_chunks(n, BLOCK_SIZE * 4, K);
    auto local_BA = std::vector<size_type>(K * num_chunks);
#pragma omp parallel for num_threads(num_threads) schedule(static, chunk_size)
    for (auto i = size_type{}; i < n; i++)
      if (check(i))
        local_BA[i / chunk_size * K + S[i]]++;
    return local_BA;
  }

  template <typename size_type>
  auto
  get_bucket(const std::ranges::random_access_range auto& S, size_type K) {
    auto n = (size_type)S.size();
    // try to compuate bucket in parallel with memory
    // less than induce_sort(4 * BLOCK_SIZE)
    auto [num_chunks, chunk_size, num_threads]
      = split_into_chunks(n, 4 * BLOCK_SIZE, K);
    auto local_BA = get_local_bucket(S, K, [](size_type) { return true; });

    auto BA = std::vector<size_type>{};
    if (num_chunks == 1) {
      BA = std::move(local_BA);
    } else {
      BA = std::vector<size_type>(K, 0);
#pragma omp parallel for num_threads(num_threads)
      for (auto j = size_type{}; j < K; j++)
        for (auto cid = size_type{}; cid < num_chunks; cid++)
          BA[j] += local_BA[cid * K + j];
    }
    std::inclusive_scan(std::begin(BA), std::end(BA), std::begin(BA),
                        std::plus<size_type>{}, size_type{1});
    return BA;
  }

  template <typename size_type>
  void
  put_lms_substr(const std::ranges::random_access_range auto& S, const auto& T,
                 const std::vector<size_type>& BA_,
                 std::ranges::random_access_range auto& SA) {
    auto n = (size_type)S.size();
#pragma omp parallel for num_threads(NUM_THREADS)
    for (auto i = size_type{}; i <= n; i++) SA[i] = EMPTY<size_type>;

    auto K = (size_type)BA_.size();
    auto [num_chunks, chunk_size, num_threads]
      = split_into_chunks(n, 4 * BLOCK_SIZE, K);
    if (num_chunks == 1) {
      auto BA = get_bucket<SUFFIX_TYPE::S_TYPE>(BA_);
      for (auto i = n - 1; ~i; i--)
        if (is_LMS(T, i))
          SA[--BA[S[i]]] = i;
    } else {
      auto local_BA
        = get_local_bucket(S, K, [&T](size_type i) { return is_LMS(T, i); });

#pragma omp parallel for num_threads(num_threads)
      for (auto chr = size_type{}; chr < K; chr++) {
        auto sum = size_type{};
        for (auto cid = num_chunks - 1; ~cid; cid--) {
          auto& cnt = local_BA[cid * K + chr];
          sum += cnt;
          cnt = sum;
        }
      }

#pragma omp parallel for num_threads(num_threads) schedule(static, chunk_size)
      for (auto i = size_type{}; i < n; i++) {
        if (is_LMS(T, i)) {
          auto chr = S[i];
          auto& cnt = local_BA[i / chunk_size * K + chr];
          SA[BA_[chr] - (cnt--)] = i;
        }
      }
    }
    SA[0] = n;
  }

  template<typename char_type, typename size_type>
  void
  prepare(const size_t L, const std::ranges::random_access_range auto& S,
          std::ranges::random_access_range auto& SA, const auto& T,
          std::ranges::random_access_range auto& RB) {
    if (L >= SA.size())
      return;
    decltype(L) R = std::min(SA.size(), L + BLOCK_SIZE);

#pragma omp parallel for num_threads(INDUCE_NUM_THREADS / 2)
    for (auto i = L; i < R; i++) {
      if (SA[i] == EMPTY<size_type> or SA[i] == 0) {
        RB[i - L] = {EMPTY<char_type>, 0};
      } else {
        auto induced_idx = SA[i] - 1;
        RB[i - L] = {S[induced_idx], T[induced_idx]};
      }
    }
  }

  template <typename size_type>
  void
  update(const size_t L, const std::ranges::random_access_range auto& WB,
         std::ranges::random_access_range auto& SA) {
    if (L >= SA.size())
      return;
    decltype(L) R = std::min(SA.size(), L + BLOCK_SIZE);

#pragma omp parallel for num_threads(INDUCE_NUM_THREADS / 2)
    for (auto i = L; i < R; i++) {
      auto& [idx, val] = WB[i - L];
      if (idx != EMPTY<size_type>) {
        SA[idx] = val;
      }
    }
  }

  template<auto induce_type, typename char_type, typename size_type>
  void
  induce_impl(const std::ranges::random_access_range auto& S, const auto& T,
              const std::ranges::input_range auto& rng, const size_type L,
              std::ranges::random_access_range auto& SA,
              std::ranges::random_access_range auto& RB,
              std::ranges::random_access_range auto& WB,
              std::ranges::random_access_range auto& BA) {
    for (size_type i : rng) {
      auto induced_idx = SA[i] - 1;

      if (SA[i] != EMPTY<size_type> and SA[i] != 0) {
        auto chr = EMPTY<char_type>;
        if (auto [c, t] = RB[i - L]; c != EMPTY<char_type>) {
          if (t == induce_type)
            chr = c;
        } else if (T[induced_idx] == induce_type) {
          chr = S[induced_idx];
        }

        if (chr == EMPTY<char_type>)
          continue;

        bool is_adjacent;
        auto pos = BA[chr];
        if constexpr (induce_type == SUFFIX_TYPE::L_TYPE) {
          BA[chr] += 1;
          is_adjacent = pos < L + (BLOCK_SIZE << 1);
        } else {
          BA[chr] -= 1;
          pos--;
          is_adjacent = pos + BLOCK_SIZE >= L;
        }

        // if pos is in adjacent block -> directly write it
        // otherwise, write it to WB
        if (is_adjacent) {
          SA[pos] = induced_idx;
          WB[i - L].first = EMPTY<size_type>;
        } else {
          WB[i - L] = {pos, induced_idx};
        }
      }
    }
  }

  template<auto induce_type, typename char_type, typename size_type>
  void
  induce(const std::ranges::random_access_range auto& S, const auto& T,
         std::vector<size_type>& BA,
         std::ranges::random_access_range auto& SA) {
    using rb_type = std::pair<char_type, uint8_t>;
    auto rb_EMPTY = std::pair{EMPTY<char_type>, 0};
    auto RBP = std::vector<rb_type>(BLOCK_SIZE, rb_EMPTY);
    auto RBI = std::vector<rb_type>(BLOCK_SIZE, rb_EMPTY);

    using wb_type = std::pair<size_type, size_type>;
    auto wb_EMPTY = std::pair{EMPTY<size_type>, EMPTY<size_type>};
    auto WBI = std::vector<wb_type>(BLOCK_SIZE, wb_EMPTY);
    auto WBU = std::vector<wb_type>(BLOCK_SIZE, wb_EMPTY);

    // views
    constexpr auto iter_view = [] {
      if constexpr (induce_type == SUFFIX_TYPE::L_TYPE) {
        return std::views::all;
      } else {
        return std::views::reverse;
      }
    }();

    auto size = (size_type)SA.size();
    auto blocks
      = std::views::iota(size_type(0), size)
        | std::views::filter([](auto i) { return i % BLOCK_SIZE == 0; });

    // prepare for first block
    if constexpr (induce_type == SUFFIX_TYPE::L_TYPE) {
      prepare<char_type, size_type>(0, S, SA, T, RBP);
    } else {
      prepare<char_type, size_type>(size / BLOCK_SIZE * BLOCK_SIZE, S, SA, T, RBP);
    }

    auto pool = biovoltron::DeprecatedThreadPool(2);
    auto stage = std::array<std::future<void>, 2>{};

    // pipeline
    for (size_type L : blocks | iter_view) {
      for (auto& s : stage)
        if (s.valid())
          s.wait();

      RBI.swap(RBP);
      WBI.swap(WBU);

      // prepare && update
      size_type P_L = L + BLOCK_SIZE;
      size_type U_L = L - BLOCK_SIZE;
      if constexpr (induce_type == SUFFIX_TYPE::S_TYPE) {
        std::swap(P_L, U_L);
      }

      stage[0] = pool.enqueue(prepare<char_type, size_type, decltype(S), decltype(SA),
                                      decltype(T), decltype(RBP)>,
                              P_L, std::ref(S), std::ref(SA), std::ref(T),
                              std::ref(RBP));
      stage[1] = pool.enqueue(update<size_type, decltype(WBU), decltype(SA)>, U_L,
                              std::ref(WBU), std::ref(SA));

      // induce
      auto rng
        = std::views::iota(L, std::min(L + BLOCK_SIZE, size)) | iter_view;
      induce_impl<induce_type, char_type>(S, T, rng, L, SA, RBI, WBI, BA);
    }
  }

  template <typename size_type>
  void
  induce_sort(const std::ranges::random_access_range auto& S, const auto& T,
              std::vector<size_type>& BA_,
              std::ranges::random_access_range auto& SA) {
    using char_type = std::remove_cvref<decltype(S.front())>::type;

    auto BA = get_bucket<SUFFIX_TYPE::L_TYPE>(BA_);
    induce<SUFFIX_TYPE::L_TYPE, char_type>(S, T, BA, SA);

    // clean SUFFIX_TYPE::S_TYPE
#pragma omp parallel for num_threads(NUM_THREADS)
    for (auto i = size_type{1}; i < SA.size(); i++)
      if (SA[i] != EMPTY<size_type> and T[SA[i]] == SUFFIX_TYPE::S_TYPE)
        SA[i] = EMPTY<size_type>;

    induce<SUFFIX_TYPE::S_TYPE, char_type>(S, T, BA_, SA);
  }

  template <typename size_type>
  auto
  is_same_substr(const std::ranges::random_access_range auto& S, const auto& T,
                 size_type x, size_type y) {
    auto n = (size_type)S.size();
    do {
      if (x == n or y == n or S[x++] != S[y++])
        return false;
    } while (not is_LMS(T, x) and not is_LMS(T, y));
    return x != n and y != n and is_LMS(T, x) and is_LMS(T, y) and S[x] == S[y];
  }

  template <typename size_type>
  auto
  name_lms_substr_left_shift(const std::ranges::random_access_range auto& S,
                             const auto& T,
                             std::ranges::random_access_range auto& SA) {
    auto n = (size_type)S.size();
    auto [num_chunks, chunk_size, num_threads]
      = split_into_chunks(n + 1, 2 * BLOCK_SIZE, 1);

    auto len = std::vector<size_type>(num_chunks, 0);
#pragma omp parallel for num_threads(num_threads) schedule(static, chunk_size)
    for (auto i = size_type{}; i <= n; i++) {
      if (is_LMS(T, SA[i])) {
        auto cid = i / chunk_size;
        auto beg = cid * chunk_size;
        std::swap(SA[beg + len[cid]++], SA[i]);
      } else
        SA[i] = EMPTY<size_type>;
    }

    auto off = std::vector<size_type>{};
    std::exclusive_scan(std::begin(len), std::end(len), std::back_inserter(off),
                        size_type{0});

#pragma omp parallel for num_threads(num_threads) schedule(dynamic, 1)
    for (auto cid = size_type{}; cid < num_chunks; cid++) {
      auto beg = cid * chunk_size;
      if (beg == off[cid])
        continue;

      for (auto i = beg; i < beg + len[cid]; i++) {
        auto j = SA[i];
        auto x = off[cid] + i - beg;

        SA[i] = EMPTY<size_type>;
        while (!__sync_bool_compare_and_swap(&SA[x], EMPTY<size_type>, j)) {}
      }
    }
    return off.back() + len.back();
  }

  template <typename size_type>
  auto
  name_lms_substr_relabel_rank(const std::ranges::random_access_range auto& S,
                               const auto& T, const auto& diff,
                               std::ranges::random_access_range auto& SA) {
    auto n1 = (size_type)diff.size();
    auto name_idx = std::vector<size_type>(NUM_THREADS);
#pragma omp parallel for num_threads(NUM_THREADS)
    for (auto i = size_type{1}; i < n1; i++) {
      auto tid = omp_get_thread_num();
      name_idx[tid] += diff[i];
    }

    exclusive_scan(std::begin(name_idx), std::end(name_idx),
                   std::begin(name_idx), 0);

#pragma omp parallel for num_threads(NUM_THREADS)
    for (auto i = size_type{1}; i < n1; i++) {
      auto tid = omp_get_thread_num();
      name_idx[tid] += diff[i];
      SA[n1 + (SA[i] >> 1)] = name_idx[tid] - 1;
    }

    return name_idx.back();
  }

  template <typename size_type>
  void
  name_lms_substr_right_shift(const std::ranges::random_access_range auto& S,
                              const auto& T, auto n1,
                              std::ranges::random_access_range auto& SA) {
    auto n = (size_type)S.size();
    auto [num_chunks, chunk_size, num_threads]
      = split_into_chunks(n - n1 + 1, 2 * BLOCK_SIZE, 1);
    auto len = std::vector<size_type>(num_chunks, 0);
#pragma omp parallel for num_threads(num_threads) schedule(static, chunk_size)
    for (auto i = n; i >= n1; i--) {
      if (SA[i] == EMPTY<size_type>)
        continue;
      auto cid = (n - i) / chunk_size;
      auto beg = cid * chunk_size;
      std::swap(SA[i], SA[n - beg - len[cid]++]);
    }

    auto off = std::vector<size_type>{};
    std::exclusive_scan(std::begin(len), std::end(len), std::back_inserter(off),
                        size_type{});

#pragma omp parallel for num_threads(num_threads) schedule(dynamic, 1)
    for (auto cid = size_type{}; cid < num_chunks; cid++) {
      auto beg = cid * chunk_size;
      if (beg == off[cid] or len[cid] == 0)
        continue;

      for (auto i = beg; i < beg + len[cid]; i++) {
        auto j = SA[n - i];
        auto x = n - (off[cid] + i - beg);

        SA[n - i] = EMPTY<size_type>;
        while (!__sync_bool_compare_and_swap(&SA[x], EMPTY<size_type>, j)) { }
      }
    }
  }

  // return the length of S1 and # kind of char in S1
  template <typename size_type>
  auto
  name_lms_substr(const std::ranges::random_access_range auto& S, const auto& T,
                  std::ranges::random_access_range auto& SA) {
    auto n = (size_type)S.size();
    auto n1 = name_lms_substr_left_shift<size_type>(S, T, SA);

    auto diff = TypeVector(n1, 0);
#pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic, BLOCK_SIZE)
    for (auto i = size_type{0}; i < n1; i++)
      diff[i] = i != 0 and not is_same_substr(S, T, SA[i - 1], SA[i]);

    auto K1 = name_lms_substr_relabel_rank<size_type>(S, T, diff, SA);

    name_lms_substr_right_shift<size_type>(S, T, n1, SA);

    // pop back tailing '\0'
    n1--;

    return std::tuple{n1, K1};
  }

  template <typename size_type>
  void
  put_lms_suffix_left_shift(const auto& T,
                            std::ranges::random_access_range auto& S1) {
    auto n = (size_type)T.size();
    auto len = std::vector<size_type>(NUM_THREADS, 0);
#pragma omp parallel for num_threads(NUM_THREADS)
    for (auto i = size_type{}; i < n; i++) {
      if (is_LMS(T, i)) {
        auto tid = omp_get_thread_num();
        len[tid]++;
      }
    }

    std::exclusive_scan(len.begin(), len.end(), len.begin(), 0);
#pragma omp parallel for num_threads(NUM_THREADS)
    for (auto i = size_type{}; i < n; i++) {
      if (is_LMS(T, i)) {
        auto tid = omp_get_thread_num();
        S1[len[tid]++] = i;
      }
    }
  }

  template <typename size_type>
  void
  put_lms_suffix_right_shift(const std::ranges::random_access_range auto& S,
                             const auto& T,
                             const std::ranges::random_access_range auto& BA_,
                             const std::ranges::random_access_range auto& S1,
                             std::ranges::random_access_range auto& SA,
                             std::ranges::random_access_range auto& SA1) {
    auto n = (size_type)S.size();
    auto n1 = (size_type)S1.size();
    auto K = (size_type)BA_.size();
    auto [num_chunks, chunk_size, num_threads]
      = split_into_chunks(n1, 4 * BLOCK_SIZE, K);
    if (num_chunks == 1) {
      auto BA = get_bucket<SUFFIX_TYPE::S_TYPE>(BA_);
      for (auto i = n1; i >= 1; i--) {
        auto j = SA1[i];
        SA1[i] = EMPTY<size_type>;
        SA[--BA[S[j]]] = j;
      }
    } else {
      auto local_BA = std::vector<size_type>(K * num_chunks);
#pragma omp parallel for num_threads(num_threads) schedule(static, chunk_size)
      for (auto i = n1; i >= 1; i--) {
        auto cid = (n1 - i) / chunk_size;
        auto chr = S[SA1[i]];
        local_BA[cid * K + chr]++;
      }

#pragma omp parallel for num_threads(num_threads)
      for (auto j = size_type{0}; j < K; j++) {
        auto sum = size_type{};
        for (auto cid = 0; cid < num_chunks; cid++) {
          auto ptr = local_BA.data() + cid * K;
          sum += ptr[j];
          ptr[j] = sum - ptr[j];
        }
      }

#pragma omp parallel for num_threads(num_threads) schedule(static, chunk_size)
      for (auto i = n1; i >= 1; i--) {
        auto cid = (n1 - i) / chunk_size;
        auto ptr = local_BA.data() + cid * K;
        auto j = SA[i];
        auto chr = S[j];
        auto x = BA_[chr] - (++ptr[chr]);
        if (j == SA[x])
          continue;

        SA[i] = EMPTY<size_type>;
        while (!__sync_bool_compare_and_swap(&SA[x], EMPTY<size_type>, j)) {}
      }
    }
  }

  template <typename size_type>
  void
  put_lms_suffix(const std::ranges::random_access_range auto& S, const auto& T,
                 const std::ranges::random_access_range auto& BA_,
                 const std::ranges::random_access_range auto& S1,
                 std::ranges::random_access_range auto& SA,
                 std::ranges::random_access_range auto& SA1) {
    put_lms_suffix_left_shift<size_type>(T, S1);

    auto n = (size_type)S.size();
    auto n1 = (size_type)S1.size();
#pragma omp parallel for num_threads(NUM_THREADS)
    for (auto i = 1; i <= n1; i++) SA1[i] = S1[SA1[i]];
    SA1[0] = n;

#pragma omp parallel for num_threads(NUM_THREADS)
    for (auto i = n1 + 1; i <= n; i++) SA[i] = EMPTY<size_type>;

    put_lms_suffix_right_shift<size_type>(S, T, BA_, S1, SA, SA1);
  }

#undef NUM_THREADS
#undef INDUCE_NUM_THREADS
} // namespace psais

} // namespace biovoltron
