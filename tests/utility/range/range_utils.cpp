#include <biovoltron/utility/range/range_utils.hpp>
#include <iostream> //debug
#include <catch.hpp>
#include <sstream>

TEST_CASE("RangeUtils::Operations - Performs various range operations", "[RangeUtils]") {
    biovoltron::RangeUtils r;
    std::ostringstream out;

    SECTION("binary_transform") {
        WARN("TODO: utility/range/range_utils.hpp requires test: binary_transform");

        std::vector<double> vec1{1.7, 2.6, 5.7, 0.0};
        std::vector<double> vec2{8.4, 9.8, 0.0, 2.7};
        std::vector<double> vec3{1.7, 9.8, 5.7, 3.8};
        std::vector<double> vec4{8.4, 2.6, 7.2, 2.7};

        auto result = r.binary_transform(vec1, vec2, [](double a, double b) { return a + b; });
        REQUIRE(result == std::vector<double>{1.7+8.4, 2.6+9.8, 5.7+0.0, 0.0+2.7});

        result = r.binary_transform(vec1, vec2, [](double a, double b) { return a - b; });
        REQUIRE(result == std::vector<double>{1.7-8.4, 2.6-9.8, 5.7-0.0, 0.0-2.7});

        result = r.binary_transform(vec1, vec2, [](double a, double b) { return a * b; });
        REQUIRE(result == std::vector<double>{1.7*8.4, 2.6*9.8, 5.7*0.0, 0.0*2.7});

        result = r.binary_transform(vec1, vec2, [](double a, double b) { return a / b; });
        REQUIRE(result == std::vector<double>{1.7/8.4, 2.6/9.8, 5.7/0.0, 0.0/2.7});

        result = r.binary_transform(vec1, vec2, [](int a, int b) { return a + b; });
        REQUIRE(result == std::vector<double>{9, 11, 5, 2});

        result = r.binary_transform(vec1, vec2, [](int a, int b) { return a - b; });
        REQUIRE(result == std::vector<double>{-7, -7, 5, -2});

        result = r.binary_transform(vec1, vec2, [](int a, int b) { return a * b; });
        REQUIRE(result == std::vector<double>{8, 18, 0, 0});

        result = r.binary_transform(vec3, vec4, [](int a, int b) { return a / b; });
        REQUIRE(result == std::vector<double>{0, 4, 0, 1});

        result = r.binary_transform(vec3, vec4, [](int a, int b) { return a % b; });
        REQUIRE(result == std::vector<double>{1, 1, 5, 1});
    }

    SECTION("index_of") {
        //std::cout << "TODO: utility/range/range_utils.hpp requires test: index_of" << std::endl;

        std::vector<int> vec{10, 20, 30, 40};
        REQUIRE(r.index_of(vec, 10) == 0);
        REQUIRE(r.index_of(vec, 30) == 2);
        REQUIRE(r.index_of(vec, 40) == 3);
        REQUIRE(r.index_of(vec, 99) == vec.size());

        std::vector<int> empty;
        REQUIRE(r.index_of(empty, 99) == 0);

        std::string s = "hello";
        REQUIRE(r.index_of(s, 'h') == 0);
        REQUIRE(r.index_of(s, 'e') == 1);
        REQUIRE(r.index_of(s, 'l') == 2);
        REQUIRE(r.index_of(s, 'x') == s.size());

        REQUIRE(r.index_of(std::vector<int>{1, 2, 3, 4}, 4) == 3);
    }

    SECTION("format_print") {
        //std::cout << "TODO: utility/range/range_utils.hpp requires test: format_print" << std::endl;

        std::vector<int> v = {1, 2, 3};

        out.str("");
        r.format_print(v, out);
        REQUIRE(out.str() == "1,2,3");

        out.str("");
        r.format_print(v, out, " | ");
        REQUIRE(out.str() == "1 | 2 | 3");

        out.str("");
        r.format_print(std::vector<int>{42}, out);
        REQUIRE(out.str() == "42");

        out.str("");
        r.format_print(std::vector<int>{}, out);
        REQUIRE(out.str().empty());

        out.str("");
        std::string s = "abc";
        r.format_print(s, out);
        REQUIRE(out.str() == "a,b,c");
    }
}

