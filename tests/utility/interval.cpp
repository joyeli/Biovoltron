#include <biovoltron/utility/interval.hpp>
#include <catch.hpp>

using namespace biovoltron;

TEST_CASE("Construction") {
  SECTION("Normal interval (begin < end)") {
    Interval intvl("chr1", 2, 10, '+');
    CHECK(intvl.chrom == "chr1");
    CHECK(intvl.begin == 2);
    CHECK(intvl.end == 10);
    CHECK(intvl.strand == '+');
    CHECK(intvl.size() == 8);
    CHECK_FALSE(intvl.empty());
  }

  GIVEN("Construct with no strand info") {
    Interval intvl("chr1", 2, 10);
    CHECK(intvl.chrom == "chr1");
    CHECK(intvl.begin == 2);
    CHECK(intvl.end == 10);
    CHECK(intvl.size() == 8);
    CHECK_FALSE(intvl.empty());
    THEN("Default is forward") {
      CHECK(intvl.strand == '+');
    }
  }

  SECTION("Interval with size 0 (begin = end)") {
    Interval intvl("chr1", 2, 2, '-');
    CHECK(intvl.chrom == "chr1");
    CHECK(intvl.begin == 2);
    CHECK(intvl.end == 2);
    CHECK(intvl.strand == '-');
    CHECK(intvl.size() == 0);
    CHECK(intvl.empty());
  }

  SECTION("Invalid interval (begin > end)") {
    CHECK_THROWS_AS(Interval("chr1", 10, 2, '+'), std::invalid_argument);
  }

  SECTION("Invalid strand symbol") {
    CHECK_THROWS_AS(Interval("chr1", 10, 2, '!'), std::invalid_argument);
  }

  SECTION("Default constuction") {
    Interval intvl;
    CHECK(intvl.chrom == "");
    CHECK(intvl.begin == 0);
    CHECK(intvl.end == 0);
    CHECK(intvl.strand == '+');
  }
}

TEST_CASE("Parse needed constuction") {
  SECTION("Normal input: forward") {
    Interval intvl("+chr1:10-2,000");
    CHECK(intvl.chrom == "chr1");
    CHECK(intvl.begin == 10);
    CHECK(intvl.end == 2000);
    CHECK(intvl.strand == '+');
  }

  SECTION("Normal input: reverse") {
    Interval intvl("-chr1:10-2,000");
    CHECK(intvl.chrom == "chr1");
    CHECK(intvl.begin == 10);
    CHECK(intvl.end == 2000);
    CHECK(intvl.strand == '-');
  }

  SECTION("Normal input: default forward") {
    Interval intvl("chr1:10-2,000");
    CHECK(intvl.chrom == "chr1");
    CHECK(intvl.begin == 10);
    CHECK(intvl.end == 2000);
    CHECK(intvl.strand == '+');
  }

  SECTION("Provide chromosome only") {
    Interval intvl("chr1");
    CHECK(intvl.chrom == "chr1");
    CHECK(intvl.begin == 0);
    CHECK(intvl.end == std::numeric_limits<std::uint32_t>::max());
    CHECK(intvl.strand == '+');
  }

  SECTION("Provide chromosome only (with strand symbol)") {
    Interval intvl("-chr1");
    CHECK(intvl.chrom == "chr1");
    CHECK(intvl.begin == 0);
    CHECK(intvl.end == std::numeric_limits<std::uint32_t>::max());
    CHECK(intvl.strand == '-');
  }

  SECTION("No end but has end of chromosome sign") {
    Interval intvl("chr1:13+");
    CHECK(intvl.chrom == "chr1");
    CHECK(intvl.begin == 13);
    CHECK(intvl.end == std::numeric_limits<std::uint32_t>::max());
    CHECK(intvl.strand == '+');
  }

  SECTION("No end but has end of chromosome sign (with strand symbol)") {
    Interval intvl("-chr1:13+");
    CHECK(intvl.chrom == "chr1");
    CHECK(intvl.begin == 13);
    CHECK(intvl.end == std::numeric_limits<std::uint32_t>::max());
    CHECK(intvl.strand == '-');
  }

  SECTION("No end and no end of chromosome sign") {
    Interval intvl("chr1:13");
    CHECK(intvl.chrom == "chr1");
    CHECK(intvl.begin == 13);
    CHECK(intvl.end == 14);
    CHECK(intvl.strand == '+');
  }

  SECTION("No end and no end of chromosome sign (with strand symbol)") {
    Interval intvl("-chr1:13");
    CHECK(intvl.chrom == "chr1");
    CHECK(intvl.begin == 13);
    CHECK(intvl.end == 14);
    CHECK(intvl.strand == '-');
  }
}

TEST_CASE("Overlap") {
  SECTION("Intervals that overlaps") {
    Interval intvl1("chr1", 100, 150);
    Interval intvl2("chr1", 120, 200);
    CHECK(intvl1.overlaps(intvl2));
  }

  SECTION("Intervals with begin >= other.end") {
    Interval intvl1("chr1", 100, 150);
    Interval intvl2("chr1", 70, 100);
    CHECK_FALSE(intvl1.overlaps(intvl2));
  }

  SECTION("Intervals with other.begin >= end") {
    Interval intvl1("chr1", 100, 150);
    Interval intvl2("chr1", 150, 200);
    CHECK_FALSE(intvl1.overlaps(intvl2));
  }

  SECTION("Intervals with different chromosome don't overlap") {
    Interval intvl1("chr1", 2, 10);
    Interval intvl2("chr2", 2, 10);
    CHECK_FALSE(intvl1.overlaps(intvl2));
  }

  SECTION("Intervals with different strand don't overlap") {
    Interval intvl1("chr1", 2, 10, '+');
    Interval intvl2("chr2", 2, 10, '-');
    CHECK_FALSE(intvl1.overlaps(intvl2));
  }
}

TEST_CASE("Contain") {
  SECTION("Interval 1 contains interval 2") {
    Interval intvl1("chr1", 100, 150);
    Interval intvl2("chr1", 120, 130);
    CHECK(intvl1.contains(intvl2));
  }

  SECTION("Begin > other.begin") {
    Interval intvl1("chr1", 100, 150);
    Interval intvl2("chr1", 70, 120);
    CHECK_FALSE(intvl1.contains(intvl2));
  }

  SECTION("End < other.end") {
    Interval intvl1("chr1", 100, 130);
    Interval intvl2("chr1", 150, 200);
    CHECK_FALSE(intvl1.contains(intvl2));
  }

  SECTION("Intervals with different chromosome doesn't contain eachother") {
    Interval intvl1("chr1", 100, 150);
    Interval intvl2("chr2", 120, 130);
    CHECK_FALSE(intvl1.contains(intvl2));
  }

  SECTION("Intervals with different stand doesn't contain eachother") {
    Interval intvl1("chr1", 100, 150, '+');
    Interval intvl2("chr2", 120, 130, '-');
    CHECK_FALSE(intvl1.contains(intvl2));
  }
}

TEST_CASE("Span") {
  SECTION("Span with interval in the same chromosome") {
    SECTION("Same strand") {
      Interval intvl1("chr1", 100, 150, '-');
      Interval intvl2("chr1", 70, 120, '-');
      auto intvl3{intvl1.span_with(intvl2)};
      CHECK(intvl3.chrom == "chr1");
      CHECK(intvl3.begin == 70);
      CHECK(intvl3.end == 150);
      CHECK(intvl3.strand == '-');
    }
    SECTION("Different strand") {
      Interval intvl1("chr1", 100, 150, '+');
      Interval intvl2("chr1", 70, 120, '-');
      CHECK_THROWS_AS(intvl1.span_with(intvl2), std::invalid_argument);
    }
  }

  SECTION("Span with interval in other chromosome") {
    Interval intvl1("chr1", 100, 150);
    Interval intvl2("chr2", 120, 130);
    CHECK_THROWS_AS(intvl1.span_with(intvl2), std::invalid_argument);
  }
}

TEST_CASE("Expand") {
  SECTION("Normal padding") {
    Interval intvl1("chr1", 100, 150, '-');
    auto intvl2{intvl1.expand_with(50)};
    CHECK(intvl2.chrom == "chr1");
    CHECK(intvl2.begin == 50);
    CHECK(intvl2.end == 200);
    CHECK(intvl2.strand == '-');
  }

  SECTION("Padding that makes begin < 0") {
    Interval intvl1("chr1", 100, 150);
    CHECK_THROWS_AS(intvl1.expand_with(120), std::invalid_argument);
  }

  SECTION("Padding that makes end overflow") {
    auto large_number = std::numeric_limits<std::uint32_t>::max() - 5;
    Interval intvl("chr1", 100, large_number);
    CHECK_THROWS_AS(intvl.expand_with(7), std::invalid_argument);
  }
}

TEST_CASE("Comparison") {
  SECTION("Equality") {
    Interval intvl1("chr1", 100, 150, '-');
    Interval intvl2("chr1", 100, 150, '-');
    CHECK(intvl1 == intvl2);
  }
  SECTION("Inequality: different chromosome") {
    Interval intvl1("chr1", 100, 150);
    Interval intvl2("achr1", 100, 150);
    CHECK(intvl1 > intvl2);  // chr1 > achr1
  }
  SECTION("Inequality: same chromosome") {
    Interval intvl1("chr1", 100, 150);
    Interval intvl2("chr1", 90, 150);
    Interval intvl3("chr1", 100, 200);
    CHECK(intvl1 > intvl2);  // 100 > 90
    CHECK(intvl1 < intvl3);  // 150 < 200
  }
  SECTION("Inequality: different strand") {
    Interval intvl1("chr1", 100, 150, '+');
    Interval intvl2("chr1", 100, 150, '-');
    CHECK(intvl1 < intvl2);  // chr1 > achr1
  }
}

TEST_CASE("Interval: to_string") {
  auto iv = Interval{"chr1", 100, 150, '-'};
  CHECK(iv.to_string() == "-chr1:100-150");
}
