#include <biovoltron/utility/read/read_filter.hpp>
#include <iostream> //debug
#include <catch.hpp>

using namespace biovoltron;

TEST_CASE("ReadFilter::Filtering - Filters reads based on various criteria", "[ReadUtils]") {
  SECTION("Mapping Quality Read Filter") {
    MappingQualityReadFilter filter;

    SamRecord<> good_read;
    good_read.mapq = 30;
    REQUIRE(filter(good_read) == false);

    SamRecord<> bad_read;
    bad_read.mapq = 10;
    REQUIRE(filter(bad_read) == true);
  }
  SECTION("Duplicate ReadFilter") {
    DuplicateReadFilter filter;

    SamRecord<> dup_read;
    dup_read.flag = SamUtil::DUPLICATE_READ;
    REQUIRE(filter(dup_read) == true);

    SamRecord<> unique_read;
    unique_read.flag = 0;
    REQUIRE(filter(unique_read) == false);
  }
  SECTION("Secondary Alignment Read Filter") {
    SecondaryAlignmentReadFilter filter;

    SamRecord<> secondary;
    secondary.flag = SamUtil::SECONDARY_ALIGNMENT;
    REQUIRE(filter(secondary) == true);

    SamRecord<> primary;
    primary.flag = 0;
    REQUIRE(filter(primary) == false);
  }
  SECTION("Minimum Length Read Filter") {
    MinimumLengthReadFilter filter;

    SamRecord<> short_read;
    short_read.seq = "ACGT"; // record.size() < 10
    REQUIRE(filter(short_read) == true);

    SamRecord<> long_read;
    long_read.seq = "ACGTACGTACGT";
    REQUIRE(filter(long_read) == false);
  }
  SECTION("Mate On Same Contig Read Filter") {
    MateOnSameContigReadFilter filter;

    SamRecord<> same_contig;
    same_contig.rnext = "=";
    REQUIRE(filter(same_contig) == false);

    SamRecord<> different_contig;
    different_contig.rnext = "chr2";
    REQUIRE(filter(different_contig) == true);
  }
}