#include <biovoltron/file_io/core/tuple.hpp>
#include <biovoltron/file_io/core/record.hpp>
#include <catch.hpp>

using namespace biovoltron;

struct TwentyOne : Record {
  int a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q,
      r, s, t, u;
};

struct Two : Record{
  char c;
  int i;
};

TEST_CASE("Tuple::to_tuple - Converts a Record to a tuple", "[FileIOcore]") {

  SECTION("Basic use") {
    auto record = Two{{}, 'a', 1};
    auto [c, i] = to_tuple(record);
    CHECK(c == 'a');
    CHECK(i == 1);
  }

  SECTION("Maximum field number") {
    auto record = TwentyOne{{}, 0, 1, 2, 3, 4, 5, 6, 7, 
      8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
    auto [a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, 
      p, q, r, s, t, u] = to_tuple(record);
    CHECK(a == 0);
    CHECK(b == 1);
    CHECK(c == 2);
    CHECK(d == 3);
    CHECK(e == 4);
    CHECK(f == 5);
    CHECK(g == 6);
    CHECK(h == 7);
    CHECK(i == 8);
    CHECK(j == 9);
    CHECK(k == 10);
    CHECK(l == 11);
    CHECK(m == 12);
    CHECK(n == 13);
    CHECK(o == 14);
    CHECK(p == 15);
    CHECK(q == 16);
    CHECK(r == 17);
    CHECK(s == 18);
    CHECK(t == 19);
    CHECK(u == 20);
  }
}
