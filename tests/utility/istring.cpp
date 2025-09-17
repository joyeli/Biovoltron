#include <biovoltron/utility/istring.hpp>
#include <bitset>
#include <catch.hpp>
#include <sstream>

using namespace biovoltron;
using namespace std::string_literals;

// We transfrom "ACGTN" to istring (integer string), 
// which was represented as 01234 internally, to allow
// a more efficient storage scheme (hashing) and ease 
// of genomic alphabet comparison.
TEST_CASE("Istring::Operations - Performs various integer string operations", "[Istring]") {
  SECTION("Regular string to istring") {
    REQUIRE(Codec::to_istring("acgt") == "0123"_s);
    REQUIRE(Codec::to_istring("ACGT") == "0123"_s);
    for (const auto& c : Codec::to_istring("bdefhijklmnopqrsuvwxyzBDEFHIJKLMNOPQRSUVWXYZ"))
      REQUIRE(c == 4);
  }

  SECTION("istring to regular string") {
    auto res = std::string{};
    std::ranges::transform(Codec::to_istring("acgt"), std::back_inserter(res), Codec::to_char);
    REQUIRE(res == "ACGT");
  }

  SECTION("A efficient storage method: hash/rhash") {
    auto dna = Codec::to_istring("aAcCgGtT");
    std::bitset<16> ans("0000010110101111");
    //                    0 0 1 1 2 2 3 3
    REQUIRE(Codec::hash(dna) == ans.to_ulong());
    REQUIRE(Codec::rhash(ans.to_ulong(), 8) == Codec::to_istring("AACCGGTT"));
  }

  SECTION("Complement") {
    std::string str = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
    std::string ans = "TNGNNNCNNNNNNNNNNNNANNNNNNTNGNNNCNNNNNNNNNNNNANNNNNN";
    for (auto& c : str) c = Codec::comp(c);
    REQUIRE(str == ans);
  }

  SECTION("Reverse complement") {
    SECTION("Regular string to regular string") {
      REQUIRE(Codec::rev_comp("atgc") == "GCAT"); 
      REQUIRE(Codec::rev_comp("ATGC") == "GCAT"); 
      REQUIRE(Codec::rev_comp("xxyy") == "NNNN"); 
    }

    SECTION("Istring to istring") {
      REQUIRE(Codec::rev_comp(Codec::to_istring("ATGC")) == Codec::to_istring("GCAT")); 
      REQUIRE(Codec::rev_comp(Codec::to_istring("xxyy")) == Codec::to_istring("NNNN")); 
    }
  }

  SECTION("operator<<") {
    std::ostringstream out;

    SECTION("Lower case") {
      out << Codec::to_istring("acgt");
      REQUIRE(out.str() == "ACGT"s);
    }

    SECTION("Upper case") {
      out << Codec::to_istring("ACGT");
      REQUIRE(out.str() == "ACGT"s);
    }

    SECTION("Other alphabet") {
      out << Codec::to_istring("bdefhijklmnopqrsuvwxyzBDEFHIJKLMNOPQRSUVWXYZ");
      for (const auto& c : out.str())
        REQUIRE(c == 'N');
    }
  }

  SECTION("operator>>") {
    auto is = istring{};
    std::istringstream in("0123");
    in >> is;
    REQUIRE(is == Codec::to_istring("NNNN"));
  }
}

TEST_CASE("Codec::Conversion - Converts between DNA characters and integers", "[Codec]") {
  SECTION("to_char") {
    REQUIRE(Codec::to_char(0u) == 'A');
    REQUIRE(Codec::to_char(1u) == 'C');
    REQUIRE(Codec::to_char(2u) == 'G');
    REQUIRE(Codec::to_char(3u) == 'T');
  }

  SECTION("to_int") {
    REQUIRE(Codec::to_int('A') == 0);
    REQUIRE(Codec::to_int('a') == 0);

    REQUIRE(Codec::to_int('C') == 1);
    REQUIRE(Codec::to_int('c') == 1);
    REQUIRE(Codec::to_int('y') == 4);
    REQUIRE(Codec::to_int('Y') == 4);
    REQUIRE(Codec::to_int('S') == 4);
    REQUIRE(Codec::to_int('s') == 4);
    REQUIRE(Codec::to_int('B') == 4);
    REQUIRE(Codec::to_int('b') == 4);

    REQUIRE(Codec::to_int('G') == 2);
    REQUIRE(Codec::to_int('g') == 2);
    REQUIRE(Codec::to_int('K') == 4);
    REQUIRE(Codec::to_int('k') == 4);

    REQUIRE(Codec::to_int('T') == 3);
    REQUIRE(Codec::to_int('t') == 3);
  }

  SECTION("is_valid") {
    std::string str = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
    for (const auto c : str) {
      switch (c) {
        case 'a': case 'A': 
        case 'c': case 'C': 
        case 'g': case 'G': 
        case 't': case 'T': REQUIRE(Codec::is_valid(c)); break;
        default: REQUIRE_FALSE(Codec::is_valid(c));
      }
    }
  }
}

TEST_CASE("Istring::Literal - Converts integer string literal", "[Istring]") {
  biovoltron::istring is = "030102030"_s;
  REQUIRE(is == biovoltron::istring{0, 3, 0, 1, 0, 2, 0, 3, 0});
}

TEST_CASE("Codec::to_int - Converts DNA character to integer", "[Codec]") {
  REQUIRE(biovoltron::Codec::to_int('A') == 0);
  REQUIRE(biovoltron::Codec::to_int('C') == 1);
  REQUIRE(biovoltron::Codec::to_int('G') == 2);
  REQUIRE(biovoltron::Codec::to_int('T') == 3);
  REQUIRE(biovoltron::Codec::to_int('N') == 4);
}

TEST_CASE("Codec::is_valid - Checks if character is a valid DNA character", "[Codec]") {
  REQUIRE(biovoltron::Codec::is_valid('A') == true);
  REQUIRE(biovoltron::Codec::is_valid('C') == true);
  REQUIRE(biovoltron::Codec::is_valid('G') == true);
  REQUIRE(biovoltron::Codec::is_valid('T') == true);
  REQUIRE(biovoltron::Codec::is_valid('N') == false);
  REQUIRE(biovoltron::Codec::is_valid('X') == false);
}

TEST_CASE("Codec::to_char - Converts integer to DNA character", "[Codec]") {
  REQUIRE(biovoltron::Codec::to_char(0) == 'A');
  REQUIRE(biovoltron::Codec::to_char(1) == 'C');
  REQUIRE(biovoltron::Codec::to_char(2) == 'G');
  REQUIRE(biovoltron::Codec::to_char(3) == 'T');
  REQUIRE(biovoltron::Codec::to_char(4) == 'N');
}

TEST_CASE("Codec::hash - Calculates hash value for DNA sequence", "[Codec]") {
  biovoltron::istring dna_iseq = "030013023"_s;
  REQUIRE(biovoltron::Codec::hash(dna_iseq) == 49611);
}

TEST_CASE("Codec::rhash - Reverses hash key to DNA sequence", "[Codec]") {
  std::size_t hash_value = 49611;
  std::size_t sequence_size = 9;
  biovoltron::istring seq = biovoltron::Codec::rhash(hash_value, sequence_size);
  REQUIRE(seq == biovoltron::istring{0, 3, 0, 0, 1, 3, 0, 2, 3});
}

TEST_CASE("Codec::rev_comp - Reverse complements a DNA sequence", "[Codec]") {
  biovoltron::istring dna_iseq = "30012303"_s;
  biovoltron::istring rev_comp_iseq = biovoltron::Codec::rev_comp(dna_iseq);
  REQUIRE(rev_comp_iseq == biovoltron::Codec::to_istring("ATACGTTA"));
}

TEST_CASE("Codec::to_string - Converts integer-based DNA sequence to string-based", "[Codec]") {
  biovoltron::istring int_seq = "21033021"_s;
  std::string string_seq = biovoltron::Codec::to_string(int_seq);
  REQUIRE(string_seq == "GCATTAGC");
}

TEST_CASE("Codec::to_istring - Converts string-based DNA sequence to integer-based", "[Codec]") {
  std::string string_seq = "TCGTAGCTGCA";
  biovoltron::istring int_iseq = biovoltron::Codec::to_istring(string_seq);
  // Replace with the actual expected integer-based sequence
  REQUIRE(int_iseq == "31230213210"_s);
}

TEST_CASE("Codec::rev_comp - Reverse complements a string-based DNA sequence", "[Codec]") {
  std::string string_seq = "TCGTCATGCTGAC";
  std::string rev_comp_iseq = biovoltron::Codec::rev_comp(string_seq);
  // Replace with the actual expected reverse complemented string-based sequence
  REQUIRE(rev_comp_iseq == "GTCAGCATGACGA"); // This is just a placeholder,
                                             // replace with the actual result
}