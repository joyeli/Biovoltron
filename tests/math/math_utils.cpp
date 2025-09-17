#include <biovoltron/math/math_utils.hpp>
#include <iostream> //debug
#include <catch.hpp>

using namespace biovoltron;

static void vector_approx_equal(const std::vector<double>& result,
  const std::vector<double>& expected, double eps = 1e-8) {
    REQUIRE(result.size() == expected.size());
    for (size_t i = 0; i < result.size(); ++i) {
        REQUIRE(result[i] == Approx(expected[i]).margin(eps));
    }
}

TEST_CASE("MathUtils::Utilities - Performs various math calculations", "[MathUtils]") {

    SECTION("get_precision") { 
      REQUIRE(MathUtils::get_precision<0>(3.5)        == Approx(4.0));
      REQUIRE(MathUtils::get_precision<2>(-3.14159)   == Approx(-3.14));
      REQUIRE(MathUtils::get_precision<3>(3.14159)    == Approx(3.142).margin(1e-4));
      REQUIRE(MathUtils::get_precision<4>(-3.14159)   == Approx(-3.1416).margin(1e-5));
      REQUIRE(MathUtils::get_precision<5>(3.1415926)  == Approx(3.14159).margin(1e-6));
      REQUIRE(MathUtils::get_precision<6>(-3.1415926) == Approx(-3.141593).margin(1e-7));
    }

    SECTION("log10_factorial") {
      for(int n = 0; n <= 12; n++){
        double expected = std::lgamma(n + 1) * std::numbers::log10e;
        REQUIRE(MathUtils::log10_factorial(n) == Approx(expected));
      }
    }

    SECTION("log10_binomial_coefficient") {
      REQUIRE(MathUtils::log10_binomial_coefficient(5, 2)  == Approx(std::log10(10)).margin(1e-12));
      REQUIRE(MathUtils::log10_binomial_coefficient(6, 3)  == Approx(std::log10(20)).margin(1e-12));
      REQUIRE(MathUtils::log10_binomial_coefficient(12, 4) == Approx(std::log10(495)).margin(1e-12));
      REQUIRE(MathUtils::log10_binomial_coefficient(5, 0)  == Approx(0));
    }

    SECTION("log_to_log10") {
      REQUIRE(MathUtils::log_to_log10(std::log(10.0))   == Approx(1.0).margin(1e-12));
      REQUIRE(MathUtils::log_to_log10(std::log(1000.0)) == Approx(3.0).margin(1e-12));
      REQUIRE(MathUtils::log_to_log10(0.0) == 0.0);
      REQUIRE(MathUtils::log_to_log10(std::log(0.01))   == Approx(-2.0).margin(1e-12));
      REQUIRE(std::isinf(MathUtils::log_to_log10(std::numeric_limits<double>::infinity())));
      REQUIRE(std::isnan(MathUtils::log_to_log10(std::numeric_limits<double>::quiet_NaN())));
    }

    SECTION("log10_gamma") {
      REQUIRE(MathUtils::log10_gamma(1)   == Approx(0.0).margin(1e-12));
      REQUIRE(MathUtils::log10_gamma(2)   == Approx(0.0).margin(1e-12));
      REQUIRE(std::isinf(MathUtils::log10_gamma(0.0))); 
      REQUIRE(std::isinf(MathUtils::log10_gamma(-1.0)));
      REQUIRE(MathUtils::log10_gamma(3)   == Approx(std::log10(2)).margin(1e-12));
      REQUIRE(MathUtils::log10_gamma(4)   == Approx(std::log10(6)).margin(1e-12));
      REQUIRE(MathUtils::log10_gamma(5)   == Approx(std::log10(24)).margin(1e-12));
      REQUIRE(MathUtils::log10_gamma(0.5) == Approx(std::log10(std::sqrt(M_PI))).margin(1e-12));
      REQUIRE(MathUtils::log10_gamma(1.5) == Approx(std::log10(0.886226925)).margin(1e-6));
      REQUIRE(MathUtils::log10_gamma(2.5) == Approx(std::log10(1.32934039)).margin(1e-6));
      REQUIRE(std::isinf(MathUtils::log10_gamma(std::numeric_limits<double>::infinity())));
      REQUIRE(std::isnan(MathUtils::log10_gamma(std::numeric_limits<double>::quiet_NaN())));
    }

    SECTION("log1mexp") {
      REQUIRE(MathUtils::log1mexp(-0.01)   == Approx(std::log(1 - std::exp(-0.01))).margin(1e-12));
      REQUIRE(MathUtils::log1mexp(-0.1)    == Approx(std::log(1 - std::exp(-0.1))).margin(1e-12));
      REQUIRE(MathUtils::log1mexp(-0.6931) == Approx(std::log(1 - std::exp(-0.6931))).margin(1e-12));
      REQUIRE(MathUtils::log1mexp(-1.0)    == Approx(std::log(1 - std::exp(-1.0))).margin(1e-12));
      REQUIRE(MathUtils::log1mexp(-1e-10)  == Approx(std::log(1 - std::exp(-1e-10))).margin(1e-12));
    }

    SECTION("log10_one_minus_pow10") {
      REQUIRE(MathUtils::log10_one_minus_pow10(-2.0)
       == Approx(std::log10(1 - std::pow(10, -2.0))).margin(1e-12));
      REQUIRE(MathUtils::log10_one_minus_pow10(-1e-10)
        == Approx(std::log10(1 - std::pow(10, -1e-10))).margin(1e-12));
      REQUIRE(MathUtils::log10_one_minus_pow10(-1e-5)
        == Approx(std::log10(1 - std::pow(10, -1e-5))).margin(1e-12));
      REQUIRE(MathUtils::log10_one_minus_pow10(-20.0)
        == Approx(std::log10(1 - std::pow(10, -20.0))).margin(1e-12));
    }

    SECTION("log10_sum_log10") {
      REQUIRE(MathUtils::log10_sum_log10(5, 5)  == Approx(5 + std::log10(2)));
      REQUIRE(MathUtils::log10_sum_log10(1, -1) == Approx(1 + std::log10(1.01)));
      REQUIRE(MathUtils::log10_sum_log10(1, 1.0000000001)
        == Approx(1.0000000001 + std::log10(2)).margin(1e-12));
      REQUIRE(MathUtils::log10_sum_log10(100, -100) == Approx(100).margin(1e-12));
      REQUIRE(MathUtils::log10_sum_log10(std::vector<double> {-100.0, -1.0, -0.5})
        == Approx(std::log10(std::pow(10.0, -100.0) + std::pow(10.0, -1.0) + std::pow(10.0, -0.5))
      ).margin(1e-8));
      REQUIRE(MathUtils::log10_sum_log10(std::vector<double>{-1.0, -1.0, -1.0})
        == Approx(std::log10(0.3)).margin(1e-8));
      REQUIRE(
        std::abs(MathUtils::log10_sum_log10(std::vector<double>{std::log10(0.3), std::log10(0.7)}))
        < 1e-12
      );
    }

    SECTION("normalize_log10") {
      REQUIRE(MathUtils::log10_sum_log10(std::vector<double>{-1.0, -1.0, -1.0})
        == Approx(std::log10(0.3)));
      REQUIRE(MathUtils::log10_sum_log10(std::vector<double>{-2.0, -1.0, -0.5})
        == Approx(std::log10(std::pow(10.0, -2.0) + std::pow(10.0, -1.0) + std::pow(10.0, -0.5)))
      );
      REQUIRE(
        std::abs(MathUtils::log10_sum_log10(std::vector<double>{std::log10(0.4), std::log10(0.6)}))
        < 1e-12
      );
      REQUIRE(MathUtils::log10_sum_log10(std::vector<double>{-3.0}) == Approx(-3.0));
    }

    SECTION("dirichlet_log10_mean_weights") {
      vector_approx_equal(
        MathUtils::dirichlet_log10_mean_weights(std::vector<double>{1.0, 2.0, 3.0}),
        {
          std::log10(1.0 / 6.0),
          std::log10(2.0 / 6.0),
          std::log10(3.0 / 6.0)
        }
      );
      vector_approx_equal(
        MathUtils::dirichlet_log10_mean_weights(std::vector<double>{1.0, 1.0, 1.0}),
        {std::vector<double>(3, std::log10(1.0 / 3.0))}
      );
      vector_approx_equal(
        MathUtils::dirichlet_log10_mean_weights(std::vector<double>{1e6, 2e6, 3e6}),
        {
          std::log10(1.0 / 6.0),
          std::log10(2.0 / 6.0),
          std::log10(3.0 / 6.0)
        }
      );
      vector_approx_equal(
        MathUtils::dirichlet_log10_mean_weights(std::vector<double>{99.0}),
        {std::log10(1.0)}
      );
    }

    SECTION("scale_log_space_array_for_numerical_stability") {
      vector_approx_equal(
        MathUtils::scale_log_space_array_for_numerical_stability(std::vector<double>{-5.0, -3.0, -1.0}),
        {-4.0, -2.0, 0.0}
      );
      vector_approx_equal(
        MathUtils::scale_log_space_array_for_numerical_stability(std::vector<double>{99.0}), {0.0}
      );
      vector_approx_equal(
        MathUtils::scale_log_space_array_for_numerical_stability(std::vector<double>{3.0, 3.0, 3.0}),
        {0.0, 0.0, 0.0}
      );
      vector_approx_equal(
        MathUtils::scale_log_space_array_for_numerical_stability(std::vector<double>{10.0, -10.0, 0.0}),
        {0.0, -20.0, -10.0}
      );
    }

    SECTION("sum_log10") {
      REQUIRE(
        MathUtils::sum_log10(std::vector<double>{std::log10(0.1), std::log10(0.1)})
        == Approx(0.2).margin(1e-12)
      );
      REQUIRE(
        MathUtils::sum_log10(std::vector<double>{std::log10(0.3), std::log10(0.7)})
        == Approx(1.0).margin(1e-12)
      );
      REQUIRE(
        MathUtils::sum_log10(std::vector<double>{std::log10(1e-10), std::log10(1.0)})
        == Approx(1.0000000001).margin(1e-12)
      );
      REQUIRE(
        MathUtils::sum_log10(std::vector<double>{std::log10(1e5), std::log10(2e5)})
        == Approx(300000.0).margin(1e-6)
      );
    }

    SECTION("normalize_from_log10_to_linear_space") {
      vector_approx_equal(
        MathUtils::normalize_from_log10_to_linear_space(
          std::vector<double>{std::log10(1.0), std::log10(1.0), std::log10(1.0)}),
        {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0}
      );
      vector_approx_equal(
        MathUtils::normalize_from_log10_to_linear_space(std::vector<double>{std::log10(42.0)}), {1.0}
      );
    }

    SECTION("approximate_log10_sum_log10") {
      REQUIRE(MathUtils::approximate_log10_sum_log10(log10(0.1), log10(0.1))
        == Approx(log10(0.2)).margin(1e-6));
      REQUIRE(MathUtils::approximate_log10_sum_log10(log10(1e-100), log10(1.0))
        == Approx(log10(1.0)).margin(1e-12));
    }
}
