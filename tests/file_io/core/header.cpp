#include <biovoltron/file_io/core/header.hpp>
#include <sstream>
#include <catch.hpp>

using namespace biovoltron;

struct FooHeader : Header {
  constexpr static auto START_SYMBOLS = std::array{"ggg", "%"};
};


TEST_CASE("Header") {
  SECTION("Basic functionality test") {
    auto content = std::string{
R"(header1
header2
header3)"
  };
    auto ss = std::stringstream{content};

    auto header = Header{};
    ss >> header;
    CHECK(header.lines.size() == 3);
    CHECK(header.lines[0] == "header1");
    CHECK(header.lines[1] == "header2");
    CHECK(header.lines[2] == "header3");

    auto ss2 = std::stringstream{};
    ss2 << header;
    CHECK(ss2.str() == content);
  }

  SECTION("How to use Header struct") {
    auto content = std::string{
R"(gggheader1
%header2
content
*content)"
    };
    auto header_content = std::string{
R"(gggheader1
%header2)"
    };
    auto ss = std::stringstream{content};

    auto header = FooHeader{};
    ss >> header;
    CHECK(header.lines.size() == 2);
    CHECK(header.lines[0] == "gggheader1");
    CHECK(header.lines[1] == "%header2");

    auto ss2 = std::stringstream{};
    ss2 << header;
    CHECK(ss2.str() == header_content);
  }
}
