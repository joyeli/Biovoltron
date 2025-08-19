#include <biovoltron/utility/variant/attribute/depth_per_allele.hpp>
#include <iostream> // debug
#include <catch.hpp>
#include <vector>

using biovoltron::DepthPerAllele;

TEST_CASE("Depth per allele", "[depth_per_allele]") {

  SECTION("get_informative_alleles") {
    // NOTE: Private method — not directly testable.
    // Keeping placeholder output.
    // std::cout << "NOTE: get_informative_alleles is private. Tested indirectly via annotate()." << std::endl;
  }

  SECTION("annotate") {
    SECTION("Returns correct support counts for informative samples") {
      std::vector<std::vector<double>> likelihoods = {
        {0.1, 0.9, 0.0},  // allele 1 (confidence = 0.8)
        {0.8, 0.1, 0.1},  // allele 0 (confidence = 0.7)
        {0.0, 0.0, 1.0}   // allele 2 (confidence = 1.0)
      };
      auto result = DepthPerAllele::annotate(likelihoods);
      REQUIRE(result == std::vector<int>({1, 1, 1}));
    }

    SECTION("Skips non-informative samples") {
      std::vector<std::vector<double>> likelihoods = {
        {0.5, 0.5, 0.0},  // confidence = 0.0 → skip
        {0.6, 0.3, 0.1},  // allele 0 (confidence = 0.3)
        {0.4, 0.3, 0.3}   // confidence = 0.1 → skip
      };
      auto result = DepthPerAllele::annotate(likelihoods);
      REQUIRE(result == std::vector<int>({1, 0, 0}));
    }

    SECTION("Handles all non-informative") {
      std::vector<std::vector<double>> likelihoods = {
        {0.5, 0.4, 0.3},
        {0.3, 0.3, 0.3}
      };
      auto result = DepthPerAllele::annotate(likelihoods);
      REQUIRE(result == std::vector<int>({0, 0, 0}));
    }

    SECTION("Threshold equality is not enough") {
      std::vector<std::vector<double>> likelihoods = {
        {0.6, 0.4, 0.0}  // confidence = 0.2 → not > threshold
      };
      auto result = DepthPerAllele::annotate(likelihoods);
      REQUIRE(result == std::vector<int>({0, 0, 0}));
    }

    SECTION("Two-allele input") {
      std::vector<std::vector<double>> likelihoods = {
        {0.9, 0.0},
        {0.1, 0.9},
        {0.51, 0.49}  // confidence = 0.02 → skip
      };
      auto result = DepthPerAllele::annotate(likelihoods);
      REQUIRE(result == std::vector<int>({1, 1}));
    }
  }
}
