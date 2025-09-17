#include <biovoltron/utility/haplotype/haplotype.hpp>
#include <iostream> //debug
#include <catch.hpp>

TEST_CASE("Haplotype::get_overlapping_events - Gets overlapping events", "[Haplotype]") {
    SECTION("get_overlapping_events") {
        biovoltron::Haplotype haplotype = {
            .event_map = {
                {20, {.location={ "chr1", 20, 21, '+' }, .ref="A", .alt="G"}},
                {30, {.location={ "chr1", 30, 31, '+' }, .ref="A", .alt="T"}},
                {40, {.location={ "chr1", 40, 41, '+' }, .ref="A", .alt="G"}},
                {10, {.location={ "chr1", 10, 11, '+' }, .ref="A", .alt="C"}},
                {50, {.location={ "chr1", 50, 51, '+' }, .ref="A", .alt="C"}}
            }
        };
        auto overlapping_events = haplotype.get_overlapping_events(30);
        REQUIRE(overlapping_events.size() == 1);
    }
}
