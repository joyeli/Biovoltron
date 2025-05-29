#pragma once
#include <biovoltron/file_io/core/record.hpp>
#include <biovoltron/utility/interval.hpp>
#include <iostream>
#include <set>
#include <unordered_set>
#include <numeric>
#include <algorithm>
#include <execution>
#include <cmath>

namespace biovoltron {

template<typename MapType>
concept is_map_type = 
  std::same_as<
    typename MapType::value_type, 
    std::pair<
      const typename MapType::key_type, 
      typename MapType::mapped_type
    >
  >;

template<class T>
concept ExpMat = is_map_type<T>; 

template <class ExpMat>
void size_normalization(ExpMat&& exp_mat) {

  auto total = 0.0;
  for (auto&& [gene, exp]: exp_mat) {
    total += exp.value;
  }

  for (auto&& [gene, exp]: exp_mat) {
    exp *= (1.0 / total);
  }

  return;
}

template <class ExpMatV>
auto fill_gene_value_and_transform(ExpMatV&& exp_mat_v) {
  std::set<std::string> gene_set;
  for (auto&& exp_mat: exp_mat_v) {
    for (auto&& [gene, exp]: exp_mat) {
      gene_set.insert(gene);
    }
  }
  
  std::vector<std::vector<double>> sample_exp_arr;

  for (auto&& exp_mat: exp_mat_v) {
    std::vector<double> tmp_exp_arr(gene_set.size());
    int idx = 0;
    for (auto&& gene_name: gene_set) {
      if (!exp_mat.contains(gene_name)) {
        exp_mat[gene_name] = typename std::remove_cvref_t<decltype(exp_mat)>::mapped_type{};
      }
      tmp_exp_arr[idx] = exp_mat[gene_name].value;
      ++idx;
    }
    sample_exp_arr.emplace_back(std::move(tmp_exp_arr));
  }

  return sample_exp_arr;
}

namespace quantile {

struct QuantileDataHelperType
{
  /// exp, idx
  std::vector<std::pair<double, int>> value_idx_pairs;
  std::vector<bool> is_ranking_drew_v;
  QuantileDataHelperType(std::vector<double>& exp_values)
	{
		for(auto idx = 0; idx < exp_values.size(); ++idx)
		{
			value_idx_pairs.emplace_back(std::make_pair(std::move(exp_values[idx]), idx));
		}
	}
};

void sort_by_rank(std::vector<QuantileDataHelperType>& q_data_v) {
  for (auto& sample: q_data_v) {
    std::sort(
      std::execution::par_unseq,
      sample.value_idx_pairs.begin(),
      sample.value_idx_pairs.end(),
      [](auto&& a, auto&& b){
        if (a.first == b.first) return a.second < b.second;
        else return a.first < b.first;
      }
    );
  }
}

void record_drew(std::vector<QuantileDataHelperType>& q_data_v) {
  /// for sorted values
  for (auto& sample: q_data_v) {
    const auto& value_pairs = sample.value_idx_pairs;
    auto& drews = sample.is_ranking_drew_v;
    for (auto idx = 0; idx < value_pairs.size(); ++idx) {
      drews.emplace_back(
        idx != value_pairs.size() && value_pairs[idx].first == value_pairs[idx + 1].first ?
        true :
        false
      );
    }
  }
}

void calc_mean_for_rank(std::vector<QuantileDataHelperType>& q_data_v) {
  auto mean = 0.0;
  auto q_data_value_size = q_data_v[0].value_idx_pairs.size();
  for (auto idx = 0; idx < q_data_value_size; ++ idx) {
    auto sum = 0.0;
    for (auto& sample: q_data_v) {
      sum += sample.value_idx_pairs[idx].first;
    }

    mean = sum / (double)(q_data_v.size());
    
    for (auto& sample: q_data_v) {
      sample.value_idx_pairs[idx].first = mean;
    }
  }
}

/// Quantile normalization 為了讓大家 distribution相同
/// 但是本來同個sample中表現量相同的annotation會因為這樣跨sample的操作變得不同
/// 所以要把他們還原成相同的數值 
/// 現在的還原順序是讓本來相同的值 變成他們的平均值
/// e.g. 本來: {10, 10, 10, 7} -> 經過轉換: {72, 15.33, 99, 12} -> 還原: {62.11, 62.11, 62.11, 12}
void replace_drew(std::vector<QuantileDataHelperType>& q_data_v) {
  auto q_data_value_size = q_data_v[0].value_idx_pairs.size();
  for (auto& sample: q_data_v) {
    auto mean = 0.0;
    const auto& drews = sample.is_ranking_drew_v;
    auto& value_pairs = sample.value_idx_pairs;
    for (auto idx = 0; idx < drews.size(); ++idx) {
      if (drews[idx] == true) {
        auto start_idx = idx;
        while (drews[idx] == true) {
          mean += value_pairs[idx].first;
          ++idx;
        }

        mean += value_pairs[idx].first;
        mean /= (double)(idx - start_idx + 1);

        for (auto j = start_idx; j <= idx; ++j) {
          value_pairs[j].first = mean;
        }
        mean = 0.0;
      }
    }
  }
}

void resume_order(std::vector<QuantileDataHelperType>& q_data_v) {
  for (auto& sample: q_data_v) {
    std::sort(
      std::execution::par_unseq,
      sample.value_idx_pairs.begin(),
      sample.value_idx_pairs.end(),
      [](auto&& a, auto&& b){
        if (a.second == b.second) return a.first < b.first;
        else return a.second < b.second;
      }
    );
  }
}

template<class ExpMatV, class QuantileNormV>
void normalize(ExpMatV&& exp_mat_v, QuantileNormV&& q_data_v) {
  /// number and order of genes in exp_mat_v and q_data_v must be the same
  for (auto sample_idx = 0; sample_idx < exp_mat_v.size(); ++sample_idx) {
    auto& exp_mat = exp_mat_v[sample_idx];
    auto& q_data = q_data_v[sample_idx];
    auto gene_idx = 0;
    for (auto&& [gene, exp]: exp_mat) {
      auto ratio = q_data.value_idx_pairs[gene_idx].first / exp.value;
      exp *= ratio;
      ++gene_idx;
    }
  }
}

template<class ExpMatV>
void quantile_normalization(ExpMatV&& exp_mat_v) {
  auto sample_exp_arr = fill_gene_value_and_transform(exp_mat_v);
  std::vector<QuantileDataHelperType> quantile_data_helper_v;
  for (auto&& exp_arr: sample_exp_arr) {
    QuantileDataHelperType tmp(exp_arr);
    quantile_data_helper_v.emplace_back(std::move(tmp));
  }
  
  sort_by_rank(quantile_data_helper_v);
  record_drew(quantile_data_helper_v);
  calc_mean_for_rank(quantile_data_helper_v);
  replace_drew(quantile_data_helper_v);
  resume_order(quantile_data_helper_v);
  normalize(exp_mat_v, quantile_data_helper_v);
}
} // namespace quantile

namespace tmm {

auto calc_library_size(const std::vector<double>& exp_arr) {
  return std::reduce(
    std::execution::par, 
    exp_arr.cbegin(),
    exp_arr.cend(),
    0.0
  );
}

auto quantile(std::vector<double> exp_arr, double p) {
  if (!std::is_sorted(exp_arr.begin(), exp_arr.end())) {
    std::sort(exp_arr.begin(), exp_arr.end());
  }
  // arr is sorted
  // interpolation='linear'
  const auto q_idx = (exp_arr.size() - 1) * p;
  const auto low_idx = std::floor(q_idx);
  const auto high_idx = std::ceil(q_idx);
  return exp_arr[low_idx] + (exp_arr[high_idx] - exp_arr[low_idx]) * (q_idx - low_idx);
}

auto pick_ref_sample(const std::vector<std::vector<double>>& sample_exp_arr) {
  auto quantiles = std::vector<double>{};
  for (auto&& exp_arr: sample_exp_arr) {
    quantiles.emplace_back(quantile(exp_arr, 0.75) / calc_library_size(exp_arr));
  }
  const auto q_mean = std::accumulate(quantiles.begin(), quantiles.end(), 0.0) / quantiles.size();
  const auto min_iter = std::min_element(quantiles.begin(), quantiles.end(), [q_mean](auto&& q1, auto&& q2){
    return std::abs(q1 - q_mean) < std::abs(q2 - q_mean);
  });
  const auto ref_idx = min_iter - quantiles.begin();
  return ref_idx;
}

auto calc_log_r(
  const std::vector<double>& obs,
  const std::vector<double>& ref,
  double n_o,
  double n_r
) {
  auto sample_size = obs.size();
  std::vector<double> log_r(sample_size);
  for (auto i = 0; i < sample_size; ++i) {
    log_r[i] = std::log2((obs[i]/n_o) / (ref[i]/n_r));
  }
  return log_r;
}

auto calc_abs_e(
  const std::vector<double>& obs,
  const std::vector<double>& ref,
  double n_o,
  double n_r
) {
  auto sample_size = obs.size();
  std::vector<double> abs_e(sample_size);
  for (auto i = 0; i < sample_size; ++i) {
    abs_e[i] = (std::log2(obs[i]/n_o) + std::log2(ref[i]/n_r)) / 2;
  }
  return abs_e;
}

auto calc_variance(
  const std::vector<double>& obs,
  const std::vector<double>& ref,
  double n_o,
  double n_r   
) {
  auto sample_size = obs.size();
  std::vector<double> var(sample_size);
  for (auto i = 0; i < sample_size; ++i) {
    var[i] = (n_o-obs[i]) / n_o / obs[i] + (n_r-ref[i]) / n_r / ref[i];
  }
  return var;
}

auto calc_f(
  const std::vector<double>& log_r,
  const std::vector<double>& variance
) {
  auto log_var_ratio_sum = 0.0;
  auto var_reciprocal_sum = 0.0;
  for (auto i = 0; i < log_r.size(); ++i) {
    log_var_ratio_sum += log_r[i] / variance[i];
    var_reciprocal_sum += 1.0 / variance[i];
  }
  return log_var_ratio_sum / var_reciprocal_sum;
}

auto calc_norm_factors_impl(
  const std::vector<double>& obs,
  const std::vector<double>& ref,
  double logratio_trim = 0.3,
  double sum_trim = 0.05,
  double A_cutoff = -1e10,
  bool do_weighting = true
) {
  auto n_r = calc_library_size(ref); 
  auto n_o = calc_library_size(obs); 
  /// M value: log ratio of expression, accounting for library size
  auto log_r = calc_log_r(obs, ref, n_o, n_r);
  /// A value: absolute expression
  auto abs_e = calc_abs_e(obs, ref, n_o, n_r);
  /// estimated asymptotic variance
  auto v = calc_variance(obs, ref, n_o, n_r);
  
  /// DEBUG
  // for (auto i: log_r) {
  //   std::cout << i << " ";
  // }
  // std::cout << "\n";
  // for (auto i: abs_e) {
  //   std::cout << i << " ";
  // }
  // std::cout << "\n";
  // for (auto i: v) {
  //   std::cout << i << " ";
  // }
  // std::cout << "\n";

  /// remove infinite and nan values, cutoff based on A
  for (int i = log_r.size() - 1; i >= 0; --i) {
    bool is_infinite_or_nan = 
      std::isinf(log_r[i]) || std::isnan(log_r[i]) ||
      std::isinf(abs_e[i]) || std::isnan(abs_e[i]) ||
      abs_e[i] < A_cutoff;

    if (is_infinite_or_nan) {
      log_r[i] = std::numeric_limits<double>::infinity();
      abs_e[i] = std::numeric_limits<double>::infinity();
      v[i] = std::numeric_limits<double>::infinity();
    }
  }
  std::erase(log_r, std::numeric_limits<double>::infinity());
  std::erase(abs_e, std::numeric_limits<double>::infinity());
  std::erase(v, std::numeric_limits<double>::infinity());

  auto abs_compare = [](auto a, auto b){
    return (std::abs(a) < std::abs(b));
  };
  if (std::abs(*max_element(abs_e.begin(), abs_e.end(), abs_compare)) < 1e-6) return 1.0;
  
  /// 把 log_r 跟 abs_e 的前後 trim 掉
  auto gene_set_size = log_r.size();
  /// 去掉前後 x% 的 gene
  /// 先計算哪些 index 會被去掉

  /// trim by M value
  auto m_low_idx = std::floor(gene_set_size * logratio_trim);
  auto m_high_idx = gene_set_size - m_low_idx;
  /// DEBUG
  // std::cout << "m low: " << m_low_idx << ", ";
  // std::cout << "m high: " << m_high_idx << "\n";

  /// trim by A value
  auto a_low_idx = std::floor(gene_set_size * sum_trim);
  auto a_high_idx = gene_set_size - a_low_idx;
  /// DEBUG
  // std::cout << "a low: " << a_low_idx << ", ";
  // std::cout << "a high: " << a_high_idx << "\n";

  /// 所以 log_r 跟 abs_e 都要先 sort完 或是知道 rank 之後才知道前後有哪些 gene

  std::set<std::size_t> keep;

  /// M value order
  /// 取 M value set 排序過後的 index
  /**
   * e.g. 
   *      unsorted idx:   0     1    2    3
   *   M value (log_r): {22.1, 5.3, 9.7, 6.8}
   *        
   *        sorted idx:   1    3    2    0  <- 我們要的 gene index (m_order)
   *      sorted array: {5.3, 6.8, 9.7, 22.1} 
   */ 
  std::vector<std::size_t> m_order(gene_set_size);
  std::iota(m_order.begin(), m_order.end(), 0);
  std::sort(
    std::execution::par_unseq,
    m_order.begin(), m_order.end(), 
    [&log_r](auto i, auto j) {
      return log_r[i] < log_r[j];
    }
  );

  /// A value order
  /// 同上
  std::vector<std::size_t> a_order(gene_set_size);
  std::iota(a_order.begin(), a_order.end(), 0);
  std::sort(
    std::execution::par_unseq,
    a_order.begin(), a_order.end(), 
    [&abs_e](auto i, auto j) {
      return abs_e[i] < abs_e[j];
    }
  );

  /// keep index that we want to retain
  /// M value keep 
  /// filter out 頭尾 x% 的 gene index
  std::set<std::size_t> keep_m;
  for (auto i = m_low_idx; i < m_high_idx; ++i) {
    keep_m.insert(m_order[i]);
  }

  /// A value keep
  /// filter out 頭尾 y% 的 gene index
  std::set<std::size_t> keep_a;
  for (auto i = a_low_idx; i < a_high_idx; ++i) {
    keep_a.insert(a_order[i]);
  }
  
  /// 取兩個 set 的交集 才是要留下來的 gene index
  std::set_intersection(
    keep_m.begin(), keep_m.end(), 
    keep_a.begin(), keep_a.end(),
    std::inserter(keep, keep.begin())
  );
  /// DEBUG
  // for (auto k: keep) {
  //   std::cout << k << " ";
  // }
  // std::cout << "\n";

  /// 把 非範圍內的 idx 去掉
  for (int i = gene_set_size - 1; i >= 0; --i) {
    if (!keep.contains(i)) {
      log_r[i] = std::numeric_limits<double>::infinity();
      v[i] = std::numeric_limits<double>::infinity();
    }
  }

  std::erase(log_r, std::numeric_limits<double>::infinity());
  std::erase(v, std::numeric_limits<double>::infinity());
  
  /// DEBUG
  // for (auto i: log_r) {
  //   std::cout << i << " ";
  // }
  // std::cout << "\n";
  // for (auto i: v) {
  //   std::cout << i << " ";
  // }
  // std::cout << "\n";

  /// if M value and A value share empty gene set, return 2^0 = 1
  if (log_r.empty() || v.empty()) return 1.0;
  
  /// calculate norm factor
  auto norm_factor = 0.0;
  if (do_weighting) {
    norm_factor = calc_f(log_r, v);
  } else {
    /// mean
    norm_factor = std::accumulate(log_r.begin(), log_r.end(), 0.0) / log_r.size();
  }

  /// DEBUG
  // std::cout << "original norm factor: " << norm_factor << std::endl;
  return std::pow(2.0, norm_factor);
}


template<class ExpMatV>
auto calc_norm_factors(ExpMatV&& exp_mat_v) {
  auto sample_size = exp_mat_v.size();
  std::vector<double> norm_factors(sample_size);
  auto sample_exp_arr = fill_gene_value_and_transform(exp_mat_v);
  auto reference_idx = pick_ref_sample(sample_exp_arr);
  for (auto i = 0; i < sample_size; ++i) {
    norm_factors[i] = calc_norm_factors_impl(
      sample_exp_arr[i], 
      sample_exp_arr[reference_idx]
    );
  }

  /// normalize norm_factors
  /// factors should multiple to one
  auto f_log_mean = 0.0;
  for (const auto& f: norm_factors) {
    f_log_mean += std::log2(f);
  }

  f_log_mean /= norm_factors.size();
  for (auto& f: norm_factors) {
    f /= std::exp2(f_log_mean);
  }

  return norm_factors;
}

template<class ExpMatV>
void normalize(ExpMatV&& exp_mat_v, const std::vector<double>& norm_factors) {
  for (auto i = 0; i < exp_mat_v.size(); ++i) {
    auto& exp_mat = exp_mat_v[i];
    auto norm_factor = norm_factors[i];
    for (auto&& kv_pair: exp_mat) {
      kv_pair.second *= (1 / norm_factor);
    }
  }
}

} // namespace tmm
}  // namespace biovoltron