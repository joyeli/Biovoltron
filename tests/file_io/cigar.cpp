#include <biovoltron/file_io/cigar.hpp>
#include <sstream>
#include <catch.hpp>

using namespace biovoltron;
using namespace std::string_literals;

// M: Alignment match (can be a sequence match or mismatch).
// S: Soft clipping (clipped sequences present in SEQ).
// D: Deletion from reference.
// I: Insertion to the reference.
// H: Hard clipping (clipped sequence NOT present in SEQ).
// N: Skipped region from the reference.
// P: Padding (silent deletion from padded reference).
// =: Sequence match.
// X: Sequence mismatch.
TEST_CASE("Cigar::Operations - Performs various CIGAR string operations", "[Cigar]") {

  SECTION("Default constructor") {
    auto cigar = Cigar{};
    CHECK(cigar == ""s);
  }

  SECTION("Constructor") {
    auto cigar = Cigar{"1M2D3I"};
    CHECK(cigar == "1M2D3I"s);
  }

  SECTION("Compact") {
    auto cigar = Cigar{"1M1M2D2D3I3I"};
    cigar.compact();
    CHECK(cigar == "2M4D6I"s);
  }

  SECTION("Emplace back") {
    auto cigar = Cigar{"1M"};
    cigar.compact();
    cigar.emplace_back(2, 'D');
    CHECK(cigar == "1M2D"s);
  }

  SECTION("Push back") {
    auto cigar = Cigar{"1M"};
    cigar.compact();
    cigar.push_back({2, 'D'});
    CHECK(cigar == "1M2D"s);
  }

  SECTION("Append") {
    auto cigar = Cigar{"1M"};
    auto other = Cigar{"2D3I"};
    cigar.append(other);
    CHECK(cigar == "1M2D3I"s);
  }

  SECTION("Swap") {
    auto cigar = Cigar{"1M2D3I"};
    auto other = Cigar{"1D"};
    cigar.swap(other);
    CHECK(cigar == "1D"s);
    CHECK(other == "1M2D3I"s);
  }

  SECTION("Reference size") {
    // Only count: M, D, N, =, X.
    auto cigar = Cigar{"1M2D3N4=5X6H"};
    CHECK(cigar.ref_size() == 15);
  }

  SECTION("Read size") {
    // Only count: M, I, S, =, X.
    auto cigar = Cigar{"1M2I3S4=5X6H"};
    CHECK(cigar.read_size() == 15);
  }

  SECTION("Begin, end, front, back, operator[]") {
    auto cigar = Cigar{"1M2D3I"};

    CHECK(static_cast<std::string>(*cigar.begin()) == "1M"s);
    CHECK_FALSE(cigar.begin() == cigar.end());

    auto front = cigar.front();
    CHECK(static_cast<std::string>(front) == "1M"s);

    auto back = cigar.back();
    CHECK(static_cast<std::string>(back) == "3I"s);

    CHECK(static_cast<std::string>(cigar[1]) == "2D"s);
  }

  SECTION("Pop front/pop back") {
    auto cigar = Cigar{"1M2D3I"};
    cigar.pop_front();
    CHECK(cigar == "2D3I"s);

    cigar.pop_back();
    CHECK(cigar == "2D"s);
  }

  SECTION("Reverse") {
    auto cigar = Cigar{"1M2D3I"};
    cigar.reverse();
    CHECK(cigar == "3I2D1M"s);
  }

  SECTION("Contains") {
    auto cigar = Cigar{"1M2D3I"};
    CHECK(cigar.contains('M'));
    CHECK_FALSE(cigar.contains('H'));
    CHECK(cigar.contains("HMS"s));
    CHECK_FALSE(cigar.contains("=NX"s));
  }

  SECTION("Size") {
    auto cigar = Cigar{"1M2D3I"};
    CHECK(cigar.size() == 3);
  }

  SECTION("Operator<<") {
    auto cigar = Cigar{"1M2D3I"};
    auto ss = std::stringstream{};
    ss << cigar;
    CHECK(ss.str() == "1M2D3I"s);
  }

  SECTION("Operator>>") {
    auto cigar = Cigar{};
    auto ss = std::stringstream{"5H5S"};
    ss >> cigar;
    CHECK(cigar == "5H5S"s);
  }

  SECTION("Unmapped read in SAMRecord") {
    auto cigar = Cigar{"*"};
    CHECK(cigar.size() == 0);
  }
}
