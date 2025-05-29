#include <biovoltron/file_io/wig.hpp>
#include <catch.hpp>
#include <filesystem>
#include <fstream>
#include <sstream>

using namespace biovoltron;
const auto data_path = std::filesystem::path{DATA_PATH};

TEST_CASE("Normal usage") {
  SECTION("Variable step") {
    std::istringstream iss{
      R"(track type=wiggle_0 name="A name" description="example file" visibility=full autoScale=off viewLimits=0.0:25.0 color=50,150,255 yLineMark=11.76 yLineOnOff=on priority=10
variableStep chrom=chr19 span=150
1000 10.1
2000 20.1
3000 30.1
)"};
    WigHeader header;
    iss >> header;
    WigVarStepRecord record;
    std::vector<WigVarStepRecord> records;
    while (iss >> record) records.emplace_back(record);

    REQUIRE(header.lines.size() == 2);
    REQUIRE(
      header.lines[0]
      == R"(track type=wiggle_0 name="A name" description="example file" visibility=full autoScale=off viewLimits=0.0:25.0 color=50,150,255 yLineMark=11.76 yLineOnOff=on priority=10)");
    REQUIRE(header.lines[1] == R"(variableStep chrom=chr19 span=150)");

    CHECK(records[0].start == 1000);
    CHECK(records[0].value == 10.1f);

    CHECK(records[1].start == 2000);
    CHECK(records[1].value == 20.1f);

    CHECK(records[2].start == 3000);
    CHECK(records[2].value == 30.1f);
  }

  SECTION("Fixed step") {
    std::istringstream iss{
      R"(track type=wiggle_0 name="A name" description="example file" visibility=full autoScale=off viewLimits=0.0:25.0 color=50,150,255 yLineMark=11.76 yLineOnOff=on priority=10
fixedStep chrom=chr19 start=1000 step=300 span=200
10.1
20.1
30.1
)"};
    WigHeader header;
    iss >> header;
    WigFixedStepRecord record;
    std::vector<WigFixedStepRecord> records;
    while (iss >> record) records.emplace_back(record);

    REQUIRE(header.lines.size() == 2);
    REQUIRE(
      header.lines[0]
      == R"(track type=wiggle_0 name="A name" description="example file" visibility=full autoScale=off viewLimits=0.0:25.0 color=50,150,255 yLineMark=11.76 yLineOnOff=on priority=10)");
    REQUIRE(header.lines[1]
            == R"(fixedStep chrom=chr19 start=1000 step=300 span=200)");

    CHECK(records[0].value == 10.1f);
    CHECK(records[1].value == 20.1f);
    CHECK(records[2].value == 30.1f);
  }
}

TEST_CASE("Read from file") {
  SECTION("Variable step") {
    std::ifstream ifs(data_path / "variableStep.wig");
    WigHeader header;
    WigVarStepRecord record;
    std::vector<WigVarStepRecord> records;

    ifs >> header;
    while (ifs >> record) records.emplace_back(record);

    REQUIRE(header.lines.size() == 9);
    REQUIRE(header.lines.back() == R"(variableStep chrom=chr19 span=150)");

    REQUIRE(records.size() == 9);
    CHECK(records[0].start == 49304701);
    CHECK(records[0].value == 10.0f);
    CHECK(records[1].start == 49304901);
    CHECK(records[1].value == 12.5f);
    CHECK(records[2].start == 49305401);
    CHECK(records[2].value == 15.0f);
    CHECK(records[3].start == 49305601);
    CHECK(records[3].value == 17.5f);
    CHECK(records[4].start == 49305901);
    CHECK(records[4].value == 20.0f);
    CHECK(records[5].start == 49306081);
    CHECK(records[5].value == 17.5f);
    CHECK(records[6].start == 49306301);
    CHECK(records[6].value == 15.0f);
    CHECK(records[7].start == 49306691);
    CHECK(records[7].value == 12.5f);
    CHECK(records[8].start == 49307871);
    CHECK(records[8].value == 10.0f);
  }

  SECTION("Fixed step") {
    std::ifstream ifs(data_path / "fixedStep.wig");
    WigHeader header;
    WigFixedStepRecord record;
    std::vector<WigFixedStepRecord> records;

    ifs >> header;
    while (ifs >> record) records.emplace_back(record);

    REQUIRE(header.lines.size() == 7);
    REQUIRE(header.lines.back()
            == R"(fixedStep chrom=chr19 start=49307401 step=300 span=200)");

    REQUIRE(records.size() == 10);
    CHECK(records[0].value == 1000.0);
    CHECK(records[1].value == 900.0);
    CHECK(records[2].value == 800.0);
    CHECK(records[3].value == 700.0);
    CHECK(records[4].value == 600.0);
    CHECK(records[5].value == 500.0);
    CHECK(records[6].value == 400.0);
    CHECK(records[7].value == 300.0);
    CHECK(records[8].value == 200.0);
    CHECK(records[9].value == 100.0);
  }
}
