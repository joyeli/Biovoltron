#include <biovoltron/file_io/core/record.hpp>
#include <sstream>
#include <catch.hpp>

using namespace biovoltron;
using namespace std::literals;

struct FooRecord : Record{
  char c;
  int i;
};

TEST_CASE("Record::IO - Reads and writes records", "[FileIOcore]") {
  SECTION("Basic use") {
    auto content = "a\t1\nb\t2\na\t1\n"s;
    auto ss = std::stringstream{content};

    auto records = std::vector<FooRecord>{};
    {
      auto record = FooRecord{};
      while (ss >> record) records.emplace_back(record);
    }

    CHECK(records[0].c == 'a');
    CHECK(records[0].i == 1);
    CHECK(records[1].c == 'b');
    CHECK(records[1].i == 2);
    CHECK(records[2].c == 'a');
    CHECK(records[2].i == 1);

    CHECK(records[0] == records[2]);

    {
      auto ss2 = std::stringstream{};
      for (const auto& record : records)
        ss2 << record;
      CHECK(ss2.str() == "a\t1\tb\t2\ta\t1\t");
    }
  }
}
