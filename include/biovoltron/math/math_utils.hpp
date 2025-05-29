#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numbers>
#include <numeric>
#include <span>

namespace biovoltron {

struct MathUtils {
  template<int N>
  static auto
  get_precision(double x) {
    constexpr auto pow10 = std::pow(10, N);
    return std::round(x * pow10) / pow10;
  }

  static auto
  log10_factorial(int n) {
    return LOG_10_FACTORIAL_CACHE[n];
  }

  static auto
  log10_binomial_coefficient(int n, int k) {
    return log10_factorial(n) - log10_factorial(k) - log10_factorial(n - k);
  }

  static auto
  log_to_log10(double ln) {
    return ln * std::numbers::log10e;
  }

  static auto
  log10_gamma(double x) -> double {
    return log_to_log10(std::lgamma(x));
  }

  static auto
  log1mexp(double a) {
    assert(a < 0);
    return (a < std::numbers::ln2) ? std::log1p(-std::exp(a)) :
                                     std::log(-std::expm1(a));
  }

  static auto
  log10_one_minus_pow10(double a) {
    assert(a < 0);
    const auto b = a * std::numbers::ln10;
    return log1mexp(b) / std::numbers::ln10;
  }

  static auto
  log10_sum_log10(double a, double b) {
    return a > b ? a + std::log10(1 + std::pow(10.0, b - a)) :
                   b + std::log10(1 + std::pow(10.0, a - b));
  }

  static auto
  log10_sum_log10(std::span<const double> log10_values) {
    const auto max_element_it = std::ranges::max_element(log10_values);
    const auto max_value = *max_element_it;
    auto sum = 1.0;
    for (auto it = log10_values.begin(); it != log10_values.end(); ++it) {
      if (it == max_element_it)
        continue;
      sum += std::pow(10.0, *it - max_value);
    }
    return max_value + (sum != 1.0 ? std::log10(sum) : 0.0);
  }

  static auto
  normalize_log10(std::span<const double> array) {
    auto result = std::vector(array.begin(), array.end());
    const auto log10_sum = log10_sum_log10(result);
    for (auto& value : result) value -= log10_sum;
    return result;
  }

  static auto
  dirichlet_log10_mean_weights(std::span<const double> alpha) {
    auto result = std::vector<double>{};
    const auto sum = std::accumulate(alpha.begin(), alpha.end(), 0.0);
    for (auto& x : alpha) result.push_back(std::log10(x / sum));
    return result;
  }

  static auto
  scale_log_space_array_for_numerical_stability(std::span<const double> array) {
    auto result = std::vector(array.begin(), array.end());
    const auto max_value = std::ranges::max(result);
    for (auto& x : result) x -= max_value;
    return result;
  }

  static auto
  sum_log10(std::span<const double> log10values) {
    return std::pow(10.0, log10_sum_log10(log10values));
  }

  static auto
  normalize_from_log10_to_linear_space(std::span<const double> array) {
    auto result = std::vector(array.begin(), array.end());
    const auto log10_sum = log10_sum_log10(result);
    for (auto& value : result) value -= log10_sum;
    for (auto& value : result) value = std::pow(10.0, value);
    return result;
  }

  static auto
  approximate_log10_sum_log10(double a, double b) -> double {
    if (a > b)
      return approximate_log10_sum_log10(b, a);
    if (const auto diff = b - a; diff < JacobianLogTable::MAX_TOLERANCE)
      return b + JacobianLogTable::get(diff);
    return b;
  }

 private:
  struct JacobianLogTable {
    constexpr static auto MAX_TOLERANCE = 8.0;
    constexpr static auto TABLE_STEP = 0.0001;
    constexpr static auto INV_STEP = 1.0 / TABLE_STEP;

    static auto
    get(double diff) -> double {
      return cache[std::round(diff * INV_STEP)];
    }

    constexpr static auto cache = [] {
      auto cache =
        std::array<double, static_cast<int>(MAX_TOLERANCE / TABLE_STEP) + 1>{};
      for (auto i = 0; i < cache.size(); i++)
        cache[i] = std::log10(1.0 + std::pow(10.0, -TABLE_STEP * i));
      return cache;
    }();
  };

  const static inline auto LOG_10_FACTORIAL_CACHE = [] {
    auto cache = std::array<double, 13>{};
    for (auto i = 0; i < cache.size(); i++) cache[i] = log10_gamma(i + 1);
    return cache;
  }();
};

}  // namespace biovoltron