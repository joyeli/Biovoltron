#include <biovoltron/utility/dna4and5.hpp>
#include <biovoltron/utility/istring.hpp>
#include <bitset>
#include <catch.hpp>
#include <sstream>
#include <iostream>

using namespace biovoltron;
using namespace std::string_literals;

TEST_CASE("Dna4::Operations - Performs various dna4 operations", "[Dna4and5]") {
  SECTION("dna4 istring") {
    REQUIRE(dna4::to_istring("acgt") == 0123_dna4);
    REQUIRE(dna4::to_istring("ACGT") == 0123_dna4);
    for (const auto& c : dna4::to_istring("bdefhijklmnopqrsuvwxyzBDEFHIJKLMNOPQRSUVWXYZ"))
      REQUIRE(c == 0);
  }
  SECTION("operation = ") {
    dna4 istring_dna4 = "ACGTacgtNzz"_dna4; // "01230123000"
    istring_dna4[0] = 'C'; // "11230123000"
    REQUIRE(istring_dna4[0]=='C');
    istring_dna4[0] = 'N'; // "01230123000"
    REQUIRE(istring_dna4[0]=='A');
  }
  SECTION("original function") {
    dna4 istring_dna4 = "ACGTacgtNzz"_dna4; // "01230123000"
    istring expected = {3, 3, 3, 0, 1, 2, 3, 0, 1, 2, 3};
    auto rev_comp_istring_dna4 = biovoltron::Codec::rev_comp(istring_dna4); // "32103210333";
    REQUIRE(rev_comp_istring_dna4==expected);
    auto hv_istring_dna4 = biovoltron::Codec::hash(istring_dna4);
    // ACGTACGTAAA
    // 0b6c6c0 = 0001101100011011000000
    REQUIRE(hv_istring_dna4==444096);
  }
  
}


TEST_CASE("Dna5::Operations - Performs various dna5 operations", "[Dna4and5]") {
  SECTION("Regular string to istring") {
    REQUIRE(dna5::to_istring("acgt") == 0123_dna5);
    REQUIRE(dna5::to_istring("ACGT") == 0123_dna5);
    for (const auto& c : dna5::to_istring("bdefhijklmnopqrsuvwxyzBDEFHIJKLMNOPQRSUVWXYZ"))
      REQUIRE(c == 4);
  }
  SECTION("operation = ") {
    dna5 istring_dna5 = "ACGTacgtNzz"_dna5; // "41230123444"
    istring_dna5[0] = 'C'; // "11230123000"
    REQUIRE(istring_dna5[0]=='C');
    istring_dna5[0] = 'N'; // "41230123000"
    REQUIRE(istring_dna5[0]=='N');
  }
  SECTION("original function") {
    dna5 istring_dna5 = "ACGTacgtNzz"_dna5; // "01230123000"
    istring expected = {4, 4, 4, 0, 1, 2, 3, 0, 1, 2, 3};
    auto rev_comp_istring_dna5 = biovoltron::Codec::rev_comp(istring_dna5);// "32103210444";
    REQUIRE(rev_comp_istring_dna5==expected);
    auto hv_istring_dna5 = biovoltron::Codec::hash(istring_dna5);
    // ACGTACGTNNN
    // 100 & 3ull is 00 for 'N'
    // so answer is 0b6c6c0
    REQUIRE(hv_istring_dna5==444096);
  }
}

