#include <biovoltron/utility/ref/hs37d5.hpp>
#include <catch.hpp>

using namespace biovoltron;

TEST_CASE("hs37d5") {
  SECTION("Get chromosome position (offset)") {
    {
      auto [name, begin, end] = Hs37d5::chr_begin_sizes[1];
      auto [chr, offset] = Hs37d5::get_chr_pos(begin+5);
      CHECK(chr == "2");
      CHECK(offset == 5);
    }
  }

  SECTION("Is valid position (not in unknown region)") {
    auto [first_begin, first_end] = Hs37d5::unknow_intervals.front();
    auto [last_begin, last_end] = Hs37d5::unknow_intervals.back();

    CHECK(Hs37d5::is_valid_pos(first_end + 5));
    CHECK_FALSE(Hs37d5::is_valid_pos(first_end - 5));
    CHECK_FALSE(Hs37d5::is_valid_pos(last_end + 10));
  }
}
