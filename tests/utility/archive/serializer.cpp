#include <biovoltron/utility/archive/serializer.hpp>
#include <catch.hpp>

using namespace biovoltron;

struct Object1 {
  int m;
  auto
  operator==(const Object1& other) const {
    return m == other.m;
  }
};

struct Object2 {
  int m;
  Object2(Object2 const&) = default;  // trivially copyable
  Object2(int x = 5) : m(x + 1) { }
  auto
  operator==(const Object2& other) const {
    return m == other.m;
  }
};

// Must be a range of trivially copyable objects.
TEMPLATE_TEST_CASE("Serializer", "", int, char, bool, float, Object1, Object2) {
  auto objs = std::vector<TestType>(5);
  {
    auto fout = std::ofstream{"temp"};
    Serializer::save(fout, objs);
  }

  auto loaded_objs = decltype(objs){};
  {
    auto fin = std::ifstream{"temp"};
    Serializer::load(fin, loaded_objs);
  }
  CHECK(objs == loaded_objs);
}
