#include <biovoltron/utility/variant/variant.hpp>
#include <iostream> //debug
#include <catch.hpp>


using namespace biovoltron;

TEST_CASE("Variant::Operations - Performs various variant operations", "[VariantUtils]") {
  SECTION("size") {
    auto variant = Variant{
      Interval{"chr1", 100, 105},
      "ATCGA", "G", {"ATCGA", "G"},
      Genotype{0, 1}, {}, 99, 40.0
    };
    REQUIRE(variant.size() == 5);
  }

  SECTION("is_snp") {
    auto snp_variant = Variant{
      Interval{"chr1", 100, 101},
      "A", "T", {"A", "T"},
      Genotype{0, 1}, {}, 99, 50.0
    };
    REQUIRE(snp_variant.is_snp());

    auto not_snp = Variant{
      Interval{"chr1", 100, 101},
      "A", "AG", {"A", "AG"},
      Genotype{0, 1}, {}, 99, 50.0
    };
    REQUIRE_FALSE(not_snp.is_snp());
  }

  SECTION("is_insertion") {
    auto insertion = Variant{
      Interval{"chr1", 100, 101},
      "A", "AG", {"A", "AG"},
      Genotype{0, 1}, {}, 99, 60.0
    };
    REQUIRE(insertion.is_insertion());

    auto not_insertion = Variant{
      Interval{"chr1", 100, 101},
      "AG", "A", {"AG", "A"},
      Genotype{0, 1}, {}, 99, 60.0
    };
    REQUIRE_FALSE(not_insertion.is_insertion());
  }

  SECTION("is_deletion") {
    auto deletion = Variant{
      Interval{"chr1", 100, 102},
      "AG", "A", {"AG", "A"},
      Genotype{0, 1}, {}, 99, 60.0
    };
    REQUIRE(deletion.is_deletion());

    auto not_deletion = Variant{
      Interval{"chr1", 100, 101},
      "A", "AG", {"A", "AG"},
      Genotype{0, 1}, {}, 99, 60.0
    };
    REQUIRE_FALSE(not_deletion.is_deletion());
  }

  SECTION("to_string") {
    auto variant = Variant{
      Interval{"chr1", 99, 100},
      "A", "G", {"A", "G"},
      Genotype{0, 1}, {0, 60, 600}, 99, 42.5
    };
    auto output = variant.to_string();

    REQUIRE_THAT(output, Catch::Matchers::StartsWith("chr1\t100\t.\tA\tG\t42.5"));
    REQUIRE(output.find("GT:GQ:PL\t0|1:99:0,60,600") != std::string::npos);
  }
}
