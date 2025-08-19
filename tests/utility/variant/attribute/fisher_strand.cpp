#include <biovoltron/utility/variant/attribute/fisher_strand.hpp>
#include <catch.hpp>
#include <cmath>  // for std::isnan

using biovoltron::FisherStrand;

TEST_CASE("Fisher strand", "[fisher_strand]") {

  SECTION("annotate returns reasonable Phred-scaled p-values") {
    double fs_score = FisherStrand::annotate(10, 10, 0, 20);  // ref: 10+10, alt: 0+20
    REQUIRE(fs_score > 20);  

    fs_score = FisherStrand::annotate(10, 10, 10, 10);
    REQUIRE(fs_score < 1e-3);  

    fs_score = FisherStrand::annotate(1, 0, 0, 5);
    REQUIRE(std::isfinite(fs_score));

  }

  SECTION("annotate downsampling logic works for large counts") {
    
    const int a = 1000, b = 1000, c = 500, d = 500;
    const double fs_score_large = FisherStrand::annotate(a, b, c, d);
    const double fs_score_scaled = FisherStrand::annotate(a / 10, b / 10, c / 10, d / 10);

    REQUIRE(std::abs(fs_score_large - fs_score_scaled) < 2.0);
  }
}
