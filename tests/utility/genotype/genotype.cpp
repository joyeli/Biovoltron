#include <biovoltron/utility/genotype/genotype.hpp>
#include <iostream> //debug
#include <catch.hpp>
#include <sstream>
TEST_CASE("Genotype::IO - Writes Genotype to stream", "[GenotypeUtils]") {

    SECTION("operator<<") {
        const auto genotype = biovoltron::Genotype{1, 2};
        std::ostringstream ss;
        ss << genotype;
        REQUIRE(ss.str() == "1|2");
    }
}
