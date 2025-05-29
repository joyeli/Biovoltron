/// tests/utility/expression/normalization.cpp
#include <biovoltron/utility/expression/normalization.hpp>
#include <biovoltron/algo/align/tailor/alignment.hpp>
#include <catch.hpp>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <concepts>

using namespace biovoltron;

struct SimExp {
  double value;

  SimExp& operator+=(const SimExp& rhs) {
    value += rhs.value;
    return *this;
  }

  SimExp& operator*=(const double val) {
    value *= val;
    return *this;
  }
};

SimExp operator+(const SimExp& lhs, const SimExp& rhs) {
  return SimExp(lhs) += rhs; 
}

SimExp operator*(const SimExp& lhs, const double val) {
  return SimExp(lhs) *= val;
}

SimExp operator*(const double val, const SimExp& rhs) {
  return SimExp(rhs) *= val;
}  

std::vector<std::map<std::string, SimExp>> samples = {
  {
    {"gene1", {15.0}}, 
    {"gene2", {12.5}}, 
    {"gene3", {22.5}}, 
    {"gene4", {1.25}}, 
    {"gene5", {8.75}} /// total 60
  },
  {
    {"gene1", {100.0}},
    {"gene3", {155.0}},
    {"gene4", {12.25}},
    {"gene6", {71.25}},
    {"gene9", {36.5}}, /// total 375
  },
  {
    {"gene1", {4.05}},
    {"gene3", {8.95}},
    {"gene4", {0.375}},
    {"gene5", {2.625}},
    {"gene8", {6.5}},
    {"gene9", {7.5}}, /// total 30 
  }
};

template<class T>
requires std::floating_point<T>
bool compare_float(T x, T y) {
  return std::abs(x - y) < 1e-9 ? true : false;
}

TEST_CASE("size normalization") {
  auto samples_copy = samples;
  auto& sample1 = samples_copy[0];
  size_normalization(sample1);
  CHECK(compare_float(sample1["gene1"].value, 0.25));
  CHECK(compare_float(sample1["gene2"].value, 12.5 / 60.0));
  CHECK(compare_float(sample1["gene3"].value, 0.375));
  CHECK(compare_float(sample1["gene4"].value, 1.25 / 60.0));
  CHECK(compare_float(sample1["gene5"].value, 8.75 / 60.0));

  auto& sample2 = samples_copy[1];
  size_normalization(sample2);
  CHECK(compare_float(sample2["gene1"].value, 100.0 / 375.0));
  CHECK(compare_float(sample2["gene3"].value, 155.0 / 375.0));
  CHECK(compare_float(sample2["gene4"].value, 12.25 / 375.0));
  CHECK(compare_float(sample2["gene6"].value, 0.19));
  CHECK(compare_float(sample2["gene9"].value, 36.5 / 375.0));

  auto& sample3 = samples_copy[2];
  size_normalization(sample3);
  CHECK(compare_float(sample3["gene1"].value, 0.135));
  CHECK(compare_float(sample3["gene3"].value, 8.95 / 30.0));
  CHECK(compare_float(sample3["gene4"].value, 0.0125));
  CHECK(compare_float(sample3["gene5"].value, 0.0875));
  CHECK(compare_float(sample3["gene8"].value, 6.5 / 30.0));
  CHECK(compare_float(sample3["gene9"].value, 0.25));
}


TEST_CASE("fill gene value and transform to sample expression array") {
  auto samples_copy = samples;
  std::vector<std::vector<double>> samples_exp_arr_ans = {
    {15.0, 12.5, 22.5, 1.25, 8.75, 0.0, 0.0, 0.0},
    {100.0, 0.0, 155.0, 12.25, 0.0, 71.25, 0.0, 36.5},
    {4.05, 0.0, 8.95, 0.375, 2.625, 0.0, 6.5, 7.5}
  };
  
  auto sample_exp_arr = fill_gene_value_and_transform(samples_copy);

  CHECK(std::equal(
    sample_exp_arr.begin(), 
    sample_exp_arr.end(), 
    samples_exp_arr_ans.begin(),
    [](auto&& v1, auto&& v2) {
      return std::equal(v1.begin(), v1.end(), v2.begin(), [](auto x, auto y){
        return compare_float(x, y);
      });
    }
  ));
  
  auto& sample1 = samples_copy[0];
  CHECK(compare_float(sample1["gene1"].value, 15.0));
  CHECK(compare_float(sample1["gene2"].value, 12.5));
  CHECK(compare_float(sample1["gene3"].value, 22.5));
  CHECK(compare_float(sample1["gene4"].value, 1.25));
  CHECK(compare_float(sample1["gene5"].value, 8.75));
  CHECK(compare_float(sample1["gene6"].value, 0.0));
  CHECK(compare_float(sample1["gene8"].value, 0.0));
  CHECK(compare_float(sample1["gene9"].value, 0.0));

  auto& sample2 = samples_copy[1];
  CHECK(compare_float(sample2["gene1"].value, 100.0));
  CHECK(compare_float(sample2["gene2"].value, 0.0));
  CHECK(compare_float(sample2["gene3"].value, 155.0));
  CHECK(compare_float(sample2["gene4"].value, 12.25));
  CHECK(compare_float(sample2["gene5"].value, 0.0));
  CHECK(compare_float(sample2["gene6"].value, 71.25));
  CHECK(compare_float(sample2["gene8"].value, 0.0));
  CHECK(compare_float(sample2["gene9"].value, 36.5));

  auto& sample3 = samples_copy[2];
  CHECK(compare_float(sample3["gene1"].value, 4.05));
  CHECK(compare_float(sample3["gene2"].value, 0.0));
  CHECK(compare_float(sample3["gene3"].value, 8.95));
  CHECK(compare_float(sample3["gene4"].value, 0.375));
  CHECK(compare_float(sample3["gene5"].value, 2.625));
  CHECK(compare_float(sample3["gene6"].value, 0.0));
  CHECK(compare_float(sample3["gene8"].value, 6.5));
  CHECK(compare_float(sample3["gene9"].value, 7.5));
}

TEST_CASE("quantile normalization") {
  std::vector<std::map<std::string, SimExp>> quantile_example = {
    {
      {"A", {5.0}}, 
      {"B", {2.0}}, 
      {"C", {3.0}}, 
      {"D", {4.0}}, 
    },
    {
      {"A", {4.0}}, 
      {"B", {1.0}}, 
      {"C", {4.0}}, 
      {"D", {2.0}},
    },
    {
      {"A", {3.0}}, 
      {"B", {4.0}}, 
      {"C", {6.0}}, 
      {"D", {8.0}},
    }
  };

  quantile::quantile_normalization(quantile_example);
  auto& sample1 = quantile_example[0];
  CHECK(compare_float(sample1["A"].value, 5.6666666667));
  CHECK(compare_float(sample1["B"].value, 2.00));
  CHECK(compare_float(sample1["C"].value, 3.00));
  CHECK(compare_float(sample1["D"].value, 4.6666666667));

  auto& sample2 = quantile_example[1];
  CHECK(compare_float(sample2["A"].value, 5.1666666667));
  CHECK(compare_float(sample2["B"].value, 2.00));
  CHECK(compare_float(sample2["C"].value, 5.1666666667));
  CHECK(compare_float(sample2["D"].value, 3.00));

  auto& sample3 = quantile_example[2];
  CHECK(compare_float(sample3["A"].value, 2.00));
  CHECK(compare_float(sample3["B"].value, 3.00));
  CHECK(compare_float(sample3["C"].value, 4.6666666667));
  CHECK(compare_float(sample3["D"].value, 5.6666666667));
}

TEST_CASE("tmm calculation utilities") {
  std::vector<std::vector<double>> samples_exp_arr = {
    {15.0, 12.5, 22.5, 1.25, 8.75, 0.0, 0.0, 0.0},
    {100.0, 0.0, 155.0, 12.25, 0.0, 71.25, 0.0, 36.5},
    {4.05, 0.0, 8.95, 0.375, 2.625, 0.0, 6.5, 7.5}
  };

  SECTION("calculation of library size") {
    auto sample1_library_size = tmm::calc_library_size(samples_exp_arr[0]);
    auto sample2_library_size = tmm::calc_library_size(samples_exp_arr[1]);
    auto sample3_library_size = tmm::calc_library_size(samples_exp_arr[2]);
    
    CHECK(compare_float(sample1_library_size, 60.0));
    CHECK(compare_float(sample2_library_size, 375.0));
    CHECK(compare_float(sample3_library_size, 30.0));
  }

  SECTION("calculation of quantile value") { 
    auto sample1_0_quantile = tmm::quantile(samples_exp_arr[0], 0.0);
    auto sample1_25_quantile = tmm::quantile(samples_exp_arr[0], 0.25);
    auto sample1_50_quantile = tmm::quantile(samples_exp_arr[0], 0.5);
    auto sample1_75_quantile = tmm::quantile(samples_exp_arr[0], 0.75);
    auto sample1_100_quantile = tmm::quantile(samples_exp_arr[0], 1);

    CHECK(compare_float(sample1_0_quantile, 0.0));
    CHECK(compare_float(sample1_25_quantile, 0.0));
    CHECK(compare_float(sample1_50_quantile, 5.0));
    CHECK(compare_float(sample1_75_quantile, 13.125));
    CHECK(compare_float(sample1_100_quantile, 22.5));
  }

  SECTION("choose the best reference sample according to the closest 0.75 quantile value to the mean") {
    auto ref_idx = tmm::pick_ref_sample(samples_exp_arr);
    CHECK(ref_idx == 0);
  }

  SECTION("calculation of M value (log_r)") {
    auto obs_sample = samples_exp_arr[1];
    auto ref_sample = samples_exp_arr[0];
    auto n_o = tmm::calc_library_size(obs_sample); 
    auto n_r = tmm::calc_library_size(ref_sample); 
    auto log_r = tmm::calc_log_r(obs_sample, ref_sample, n_o, n_r);
    
    CHECK(compare_float(log_r[0], 0.0931094044));
    CHECK(std::isinf(log_r[1]));
    CHECK(compare_float(log_r[2], 0.1404151192));
    CHECK(compare_float(log_r[3], 0.6489255595));
    CHECK(std::isinf(log_r[4]));
    CHECK(std::isinf(log_r[5]));
    CHECK(std::isnan(log_r[6]));
    CHECK(std::isinf(log_r[7]));
  }

  SECTION("calculation of A value (abs_e)") {
    auto obs_sample = samples_exp_arr[1];
    auto ref_sample = samples_exp_arr[0];
    auto n_o = tmm::calc_library_size(obs_sample); 
    auto n_r = tmm::calc_library_size(ref_sample); 
    auto abs_e = tmm::calc_abs_e(obs_sample, ref_sample, n_o, n_r);
    
    CHECK(compare_float(abs_e[0], -1.9534452978));
    CHECK(std::isinf(abs_e[1]));
    CHECK(compare_float(abs_e[2], -1.3448299397));
    CHECK(compare_float(abs_e[3], -5.260499721));
    CHECK(std::isinf(abs_e[4]));
    CHECK(std::isinf(abs_e[5]));
    CHECK(std::isinf(abs_e[6]));
    CHECK(std::isinf(abs_e[7]));
  }

  SECTION("calculation of A value (abs_e)") {
    auto obs_sample = samples_exp_arr[1];
    auto ref_sample = samples_exp_arr[0];
    auto n_o = tmm::calc_library_size(obs_sample); 
    auto n_r = tmm::calc_library_size(ref_sample); 
    auto v = tmm::calc_variance(obs_sample, ref_sample, n_o, n_r);
    
    CHECK(compare_float(v[0], 0.0573333333));
    CHECK(std::isinf(v[1]));
    CHECK(compare_float(v[2], 0.031562724));
    CHECK(compare_float(v[3], 0.8622993197));
    CHECK(std::isinf(v[4]));
    CHECK(std::isinf(v[5]));
    CHECK(std::isinf(v[6]));
    CHECK(std::isinf(v[7]));
  }

  SECTION("calculation of norm factor") {
    auto obs_sample_1 = samples_exp_arr[0];
    auto obs_sample_2 = samples_exp_arr[1];
    auto obs_sample_3 = samples_exp_arr[2];
    auto ref_sample = samples_exp_arr[0];
    auto norm_factor_1 = tmm::calc_norm_factors_impl(obs_sample_1, ref_sample);
    auto norm_factor_2 = tmm::calc_norm_factors_impl(obs_sample_2, ref_sample);
    auto norm_factor_3 = tmm::calc_norm_factors_impl(obs_sample_3, ref_sample);
    CHECK(compare_float(norm_factor_1, 1.0));
    CHECK(compare_float(norm_factor_2, 1.0986516715));
    CHECK(compare_float(norm_factor_3, 0.6));
    auto norm_factors = std::vector<double>{norm_factor_1, norm_factor_2, norm_factor_3};
    auto f_log_mean = 0.0;

    /// normalize norm_factors
    for (const auto& f: norm_factors) {
      f_log_mean += std::log2(f);
    }
  
    f_log_mean /= norm_factors.size();
    for (auto& f: norm_factors) {
      f /= std::exp2(f_log_mean);
    }

    for (auto f: norm_factors) {
      std::cout << f << " ";
    }
    std::cout << "\n";
  }

  SECTION("calculation and normalization of norm factors") {
    auto samples_copy = samples;
    auto norm_factors = tmm::calc_norm_factors(samples_copy);
    CHECK(compare_float(norm_factors[0], 1.14902526));
    CHECK(compare_float(norm_factors[1], 1.2623785225));
    CHECK(compare_float(norm_factors[2], 0.689415156));
  }

  SECTION("tmm normalization") {
    auto samples_copy = samples;
    auto norm_factors = tmm::calc_norm_factors(samples_copy);
    tmm::normalize(samples_copy, norm_factors);
    
    auto& sample1 = samples_copy[0];
    CHECK(compare_float(sample1["gene1"].value, 13.0545432912));
    CHECK(compare_float(sample1["gene2"].value, 10.878786076));
    CHECK(compare_float(sample1["gene3"].value, 19.5818149368));
    CHECK(compare_float(sample1["gene4"].value, 1.0878786076));
    CHECK(compare_float(sample1["gene5"].value, 7.6151502532));
    CHECK(compare_float(sample1["gene6"].value, 0.0));
    CHECK(compare_float(sample1["gene8"].value, 0.0));
    CHECK(compare_float(sample1["gene9"].value, 0.0));

    auto& sample2 = samples_copy[1];
    CHECK(compare_float(sample2["gene1"].value, 79.2155428936));
    CHECK(compare_float(sample2["gene2"].value, 0.0));
    CHECK(compare_float(sample2["gene3"].value, 122.784091485));
    CHECK(compare_float(sample2["gene4"].value, 9.7039040045));
    CHECK(compare_float(sample2["gene5"].value, 0.0));
    CHECK(compare_float(sample2["gene6"].value, 56.4410743117));
    CHECK(compare_float(sample2["gene8"].value, 0.0));
    CHECK(compare_float(sample2["gene9"].value, 28.9136731562));

    auto& sample3 = samples_copy[2];
    CHECK(compare_float(sample3["gene1"].value, 5.8745444811));
    CHECK(compare_float(sample3["gene2"].value, 0.0));
    CHECK(compare_float(sample3["gene3"].value, 12.9820180507));
    CHECK(compare_float(sample3["gene4"].value, 0.5439393038));
    CHECK(compare_float(sample3["gene5"].value, 3.8075751266));
    CHECK(compare_float(sample3["gene6"].value, 0.0));
    CHECK(compare_float(sample3["gene8"].value, 9.4282812659));
    CHECK(compare_float(sample3["gene9"].value, 10.878786076));
  }

  SECTION("calculation of norm factor: real case") {
    /// case from: https://www.tspweb.com/key/tmm%20%E6%A0%87%E5%87%86.html
    std::vector<double> c1(50, 10.0); /// ref_sample
    std::vector<double> c2(50, 11.0);
    std::vector<double> p1(25, 20.0);
    std::vector<double> p2(25, 21.0);
    p1.resize(50);
    p2.resize(50);
    auto norm_factor_1 = tmm::calc_norm_factors_impl(c1, c1);
    auto norm_factor_2 = tmm::calc_norm_factors_impl(c2, c1);
    auto norm_factor_3 = tmm::calc_norm_factors_impl(p1, c1);
    auto norm_factor_4 = tmm::calc_norm_factors_impl(p2, c1);

    CHECK(norm_factor_1 == 1.0);
    CHECK(norm_factor_2 == 1.0);
    CHECK(norm_factor_3 == 2.0);
    CHECK(norm_factor_4 == 2.0);

    auto norm_factors = std::vector<double>{norm_factor_1, norm_factor_2, norm_factor_3, norm_factor_4};
    auto f_log_mean = 0.0;
    for (const auto& f: norm_factors) {
      f_log_mean += std::log2(f);
    }
  
    f_log_mean /= norm_factors.size();
    for (auto& f: norm_factors) {
      f /= std::exp2(f_log_mean);
    }

    CHECK(compare_float(norm_factors[0], 0.7071067812));
    CHECK(compare_float(norm_factors[1], 0.7071067812));
    CHECK(compare_float(norm_factors[2], 1.4142135624));
    CHECK(compare_float(norm_factors[3], 1.4142135624));
  }
}
