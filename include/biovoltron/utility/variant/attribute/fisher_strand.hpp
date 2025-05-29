#pragma once

#include <biovoltron/utility/read/quality_utils.hpp>
#include <boost/math/distributions/hypergeometric.hpp>

namespace biovoltron {

class FisherStrand {
  constexpr static auto TARGET_TABLE_SIZE = 200.0;

  static auto
  fisher_test(int a, int b, int c, int d) {
    const auto N = a + b + c + d;
    const auto r = a + c;
    const auto n = c + d;
    const auto max_k = std::min(r, n);
    const auto min_k = std::max(0, r + n - N);
    const auto dist = boost::math::hypergeometric_distribution(r, n, N);
    const auto cutoff = boost::math::pdf(dist, c);
    auto prob = 0.0;
    for (auto k = min_k; k < max_k + 1; k++)
      if (const auto p = boost::math::pdf(dist, k); p <= cutoff)
        prob += p;
    return prob;
  }

 public:
  static auto
  annotate(int a, int b, int c, int d) {
    if (const auto N = a + b + c + d; N > TARGET_TABLE_SIZE * 2) {
      const auto k = N / TARGET_TABLE_SIZE;
      a /= k;
      b /= k;
      c /= k;
      d /= k;
    }
    return QualityUtils::phred_scale_error_rate(fisher_test(a, b, c, d));
  }
};

}  // namespace biovoltron