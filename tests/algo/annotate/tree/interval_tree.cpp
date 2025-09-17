#include <biovoltron/algo/annotate/tree/interval_tree.hpp>
#include <biovoltron/utility/interval.hpp>
#include <catch.hpp>

using namespace biovoltron;

TEST_CASE("IntervalTree::find - Finds overlapping intervals", "[IntervalTree]") {
  SECTION("Build interval tree and retrieve overlaping data") {
    IntervalTree<std::string> itree;
    itree.insert(5, 10, "data1");
    itree.insert(2, 13, "data2");
    itree.insert(20, 30, "data3");
    itree.index();
    auto results = itree.find(6, 9);
    // result is a vector of object, sorted by begin position
    REQUIRE(results.size() == 2);
    REQUIRE(results[0] == "data2");
    REQUIRE(results[1] == "data1");

    results = itree.find(30, 40);
    REQUIRE(results.empty());
  }

  SECTION("Insert object that has location info") {
    IntervalTree<std::string> itree;
    itree.insert(5, 10, Interval{"data1", 5, 10}.to_string());
    itree.insert(2, 13, Interval{"data2", 2, 13}.to_string());
    itree.insert(20, 30, Interval{"data3", 20, 30}.to_string());
    itree.index();
    auto results = itree.find(6, 9);
    REQUIRE(results.size() == 2);
    REQUIRE(results[0] == Interval{"data2", 2, 13});
    REQUIRE(results[1] == Interval{"data1", 5, 10});

    results = itree.find(30, 40);
    REQUIRE(results.empty());
  }

  SECTION("Throw when find() before index()") {
    IntervalTree<std::string> itree;
    itree.insert(2, 10, "data");
    REQUIRE_THROWS(itree.find(5, 9));
  }

  SECTION("A test that run through all the logic in find()") {
    IntervalTree<std::string> itree;
    // root
    itree.insert(150, 160, "data2");
    // left child
    itree.insert(32, 80, "data1");
    // left grand
    itree.insert(30, 200, "data2");
    itree.insert(29, 70, "data1");
    itree.insert(28, 70, "data1");
    itree.insert(27, 70, "data1");
    itree.insert(26, 70, "data1");
    itree.insert(25, 70, "data1");
    itree.insert(24, 70, "data1");
    itree.insert(23, 70, "data1");
    itree.insert(22, 70, "data1");
    itree.insert(21, 70, "data1");
    itree.insert(20, 70, "data1");
    itree.insert(19, 70, "data1");
    itree.insert(18, 70, "data1");
    itree.insert(17, 70, "data1");
    itree.insert(16, 70, "data1");
    itree.insert(15, 70, "data1");
    itree.insert(14, 70, "data1");
    itree.insert(13, 70, "data1");
    itree.insert(12, 70, "data1");
    itree.insert(11, 70, "data1");
    itree.insert(10, 70, "data1");
    itree.insert(9, 70, "data1");
    itree.insert(8, 70, "data1");
    itree.insert(7, 70, "data1");
    itree.insert(6, 70, "data1");
    itree.insert(5, 70, "data1");
    itree.insert(4, 70, "data1");
    itree.insert(3, 70, "data1");
    itree.insert(2, 70, "data1");
    itree.insert(1, 70, "data1");
    itree.insert(0, 70, "data1");
    // right grand
    itree.insert(68, 80, "data1");
    itree.insert(67, 80, "data1");
    itree.insert(66, 80, "data1");
    itree.insert(65, 80, "data1");
    itree.insert(64, 80, "data1");
    itree.insert(63, 80, "data1");
    itree.insert(62, 80, "data1");
    itree.insert(61, 80, "data1");
    itree.insert(60, 80, "data1");
    itree.insert(59, 80, "data1");
    itree.insert(58, 80, "data1");
    itree.insert(57, 80, "data1");
    itree.insert(56, 80, "data1");
    itree.insert(55, 80, "data1");
    itree.insert(54, 80, "data1");
    itree.insert(53, 80, "data1");
    itree.insert(52, 80, "data1");
    itree.insert(51, 80, "data1");
    itree.insert(50, 80, "data1");
    itree.insert(49, 80, "data1");
    itree.insert(48, 80, "data1");
    itree.insert(47, 80, "data1");
    itree.insert(46, 80, "data1");
    itree.insert(45, 80, "data1");
    itree.insert(44, 80, "data1");
    itree.insert(43, 80, "data1");
    itree.insert(42, 80, "data1");
    itree.insert(41, 80, "data1");
    itree.insert(40, 80, "data1");
    itree.insert(39, 80, "data1");
    itree.insert(38, 80, "data1");
    // right right subtree
    itree.insert(170, 300, "data2");
    itree.insert(168, 250, "data2");
    itree.insert(172, 250, "data2");
    itree.insert(166, 250, "data2");
    itree.insert(173, 250, "data2");
    itree.insert(165, 250, "data2");
    itree.insert(174, 250, "data2");
    itree.insert(164, 250, "data2");
    itree.insert(175, 250, "data2");
    itree.insert(163, 250, "data2");
    itree.insert(176, 250, "data2");
    itree.insert(162, 250, "data2");
    itree.insert(177, 250, "data2");
    itree.insert(161, 250, "data2");
    itree.insert(210, 250, "data1");

    itree.index();
    auto results = itree.find(100, 200);

    // all the overlaps are labeled "data2", there should be 16 of them
    REQUIRE(results.size() == 16);
    for (const auto& result : results) REQUIRE(result == "data2");
  }
}
