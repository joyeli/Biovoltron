#include <biovoltron/utility/read/quality_utils.hpp>
#include <iostream> //debug
#include <catch.hpp>

TEST_CASE("Quality utils") {

  using biovoltron::QualityUtils;

  SECTION("qual_to_error_prob") {

    REQUIRE(QualityUtils::qual_to_error_prob('!' - QualityUtils::ASCII_OFFSET) == Approx(1.0)); // Q=0
    REQUIRE(QualityUtils::qual_to_error_prob('+' - QualityUtils::ASCII_OFFSET) == Approx(0.1).epsilon(1e-10));     // Q=10
    REQUIRE(QualityUtils::qual_to_error_prob('5' - QualityUtils::ASCII_OFFSET) == Approx(0.01).epsilon(1e-10));    // Q=20
    REQUIRE(QualityUtils::qual_to_error_prob('?' - QualityUtils::ASCII_OFFSET) == Approx(0.001).epsilon(1e-10));   // Q=30
    REQUIRE(QualityUtils::qual_to_error_prob('I' - QualityUtils::ASCII_OFFSET) == Approx(0.0001).epsilon(1e-10));  // Q=40
    REQUIRE(QualityUtils::qual_to_error_prob('~' - QualityUtils::ASCII_OFFSET) == Approx(std::pow(10, ('~' - QualityUtils::ASCII_OFFSET) / -10.0)).epsilon(1e-10));
    
    for (int q = 1; q < 127; ++q) {
      auto prev = QualityUtils::qual_to_error_prob(q - 1);
      auto curr = QualityUtils::qual_to_error_prob(q);
      REQUIRE(curr <= prev); 
    }
  }

  SECTION("qual_to_error_prob_log10") {
    //  Test edge cases
    REQUIRE(QualityUtils::qual_to_error_prob_log10(0) == Approx(0.0));

    REQUIRE(QualityUtils::qual_to_error_prob_log10(10) == Approx(-1.0).epsilon(1e-10));
    REQUIRE(QualityUtils::qual_to_error_prob_log10(20) == Approx(-2.0).epsilon(1e-10));
    REQUIRE(QualityUtils::qual_to_error_prob_log10(30) == Approx(-3.0).epsilon(1e-10));
    REQUIRE(QualityUtils::qual_to_error_prob_log10(40) == Approx(-4.0).epsilon(1e-10));

    //  Test for values outside the expected range
    REQUIRE(QualityUtils::qual_to_error_prob_log10(-10) == Approx(1.0));
  }

  SECTION("qual_to_prob_log10") {
      auto val = QualityUtils::qual_to_prob_log10(0);

      REQUIRE(std::isinf(val));
      REQUIRE(val < 0);  // -inf

      REQUIRE(QualityUtils::qual_to_prob_log10(10) == Approx(std::log10(0.9)).epsilon(1e-10));
      REQUIRE(QualityUtils::qual_to_prob_log10(20) == Approx(std::log10(0.99)).epsilon(1e-10));
      REQUIRE(QualityUtils::qual_to_prob_log10(30) == Approx(std::log10(0.999)).epsilon(1e-10));
      REQUIRE(QualityUtils::qual_to_prob_log10(40) == Approx(std::log10(0.9999)).epsilon(1e-10));

  }

  SECTION("phred_scale_error_rate") {
      REQUIRE(QualityUtils::phred_scale_error_rate(1.0) == Approx(0.0).epsilon(1e-10));
      REQUIRE(QualityUtils::phred_scale_error_rate(0.1) == Approx(10.0).epsilon(1e-10));
      REQUIRE(QualityUtils::phred_scale_error_rate(0.01) == Approx(20.0).epsilon(1e-10));
      REQUIRE(QualityUtils::phred_scale_error_rate(0.001) == Approx(30.0).epsilon(1e-10));
      REQUIRE(QualityUtils::phred_scale_error_rate(0.0001) == Approx(40.0).epsilon(1e-10));

      REQUIRE(QualityUtils::phred_scale_error_rate(0.005) >
              QualityUtils::phred_scale_error_rate(0.05));
  }
}

