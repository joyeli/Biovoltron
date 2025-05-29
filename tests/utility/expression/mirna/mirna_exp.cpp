/// tests/utility/expression/mirna/mirna_exp.cpp
#include <biovoltron/utility/expression/mirna/mirna_exp.hpp>
#include <biovoltron/algo/align/tailor/alignment.hpp>
#include <catch.hpp>
#include <filesystem>
#include <sstream>
#include <iostream>
#include <algorithm>

/// mirna expression type constructor

using namespace biovoltron;
///                           value, lens
auto total_value = 10.0;
auto a_tail = mirna::TailExp { 1.0, { {18, {0.5}}, {19, {0.5}} } };
auto c_tail = mirna::TailExp { 1.0, { {18, {0.25}}, {19, {0.25}}, {20, {0.5}} } };
auto g_tail = mirna::TailExp { 1.0, { {19, {0.5}}, {20, {0.5}} } };
auto t_tail = mirna::TailExp { 1.0, { {20, {0.25}}, {21, {0.25}}, {23, {0.5}} } };
auto o_tail = mirna::TailExp { 1.0, { {22, {1.0}} } };
auto gm_exp = mirna::TailExp { 5.0, { {18, {2.5}}, {23, {1.0}}, {25, {1.5}} } };

TEST_CASE("MirExp constructor") {
  
  SECTION("empty object") {
    auto empty_mir = mirna::MirExp{};
    CHECK(empty_mir.value == 0);

    CHECK(empty_mir.tails[0].value == 0);
    CHECK(empty_mir.tails[0].lens.empty());
    CHECK(empty_mir.tails[0].lens[21].value == 0);

    CHECK(empty_mir.tails[1].value == 0);
    CHECK(empty_mir.tails[1].lens.empty());
    CHECK(empty_mir.tails[1].lens[18].value == 0);

    CHECK(empty_mir.tails[2].value == 0);
    CHECK(empty_mir.tails[2].lens.empty());
    CHECK(empty_mir.tails[2].lens[19].value == 0);

    CHECK(empty_mir.tails[3].value == 0);
    CHECK(empty_mir.tails[3].lens.empty());
    CHECK(empty_mir.tails[3].lens[22].value == 0);

    CHECK(empty_mir.tails[4].value == 0);
    CHECK(empty_mir.tails[4].lens.empty());
    CHECK(empty_mir.tails[4].lens[20].value == 0);

    CHECK(empty_mir.tails[5].value == 0);
    CHECK(empty_mir.tails[5].lens.empty());
    CHECK(empty_mir.tails[5].lens[2].value == 0);
  }
  
  SECTION("initialization list") {
    auto mir_1 = mirna::MirExp {
      total_value, // value 
      {{
        a_tail,
        c_tail,
        g_tail,
        t_tail,
        o_tail,
        gm_exp
      }}
    };
    CHECK(mir_1.value == 10.0);

    // A tail
    CHECK(mir_1.tails[0].value == 1.0);
    REQUIRE(!mir_1.tails[0].lens.empty());
    CHECK(mir_1.tails[0].lens[18].value == 0.5);
    CHECK(mir_1.tails[0].lens[19].value == 0.5);

    // C tail
    CHECK(mir_1.tails[1].value == 1.0);
    REQUIRE(!mir_1.tails[1].lens.empty());
    CHECK(mir_1.tails[1].lens[18].value == 0.25);
    CHECK(mir_1.tails[1].lens[19].value == 0.25);
    CHECK(mir_1.tails[1].lens[20].value == 0.5);
    
    // G tail
    CHECK(mir_1.tails[2].value == 1.0);
    REQUIRE(!mir_1.tails[2].lens.empty());
    CHECK(mir_1.tails[2].lens[19].value == 0.5);
    CHECK(mir_1.tails[2].lens[20].value == 0.5);

    // T tail
    CHECK(mir_1.tails[3].value == 1.0);
    REQUIRE(!mir_1.tails[3].lens.empty());
    CHECK(mir_1.tails[3].lens[20].value == 0.25);
    CHECK(mir_1.tails[3].lens[21].value == 0.25);
    CHECK(mir_1.tails[3].lens[23].value == 0.5);
    
    // Other tail
    CHECK(mir_1.tails[4].value == 1.0);
    REQUIRE(!mir_1.tails[4].lens.empty());  
    CHECK(mir_1.tails[4].lens[22].value == 1.0);

    // Genome mathcing
    CHECK(mir_1.tails[5].value == 5.0);
    REQUIRE(!mir_1.tails[5].lens.empty());  
    CHECK(mir_1.tails[5].lens[18].value == 2.5);
    CHECK(mir_1.tails[5].lens[23].value == 1.0);
    CHECK(mir_1.tails[5].lens[25].value == 1.5);
  }
  SECTION("construct from alignment") {
    auto hit1 = Hit{
      {{4, 'T'}, {1, 'C'}},
      {},
      {"chr1", 0, 10}
    };
    auto hit2 = Hit{
      {{4, 'T'}, {1, 'C'}},
      {},
      {"chr2", 10, 20}
    };
    auto hit3 = Hit{
      {{4, 'A'}, {1, 'G'}},
      {},
      {"chr3", 20, 30}
    };
    auto aln1 = Alignment{
      "seq1",
      "AACCGGTTGG",
      "!!!!!!!!!!",
      true, // forward?
      -1, // tail position
      {hit1, hit2, hit3},
      3
    };
    auto aln2 = Alignment{
      "seq2",
      "AACCGGTTGG",
      "!!!!!!!!!!",
      true, // forward?
      8, // tail position
      {hit1, hit2},
      2
    };
    auto mir_1 = mirna::MirExp::init_from_alignment(aln1);
    CHECK(mir_1.value == (1.0 / 3.0));
    CHECK(mir_1.tails[5].lens[10].value == (1.0 / 3.0)); // genome matching, len 10
    CHECK(mir_1.tails[5].lens[12].value == 0); // genome matching, len 12
    CHECK(mir_1.tails[1].lens[10].value == 0); // C tail, len 10

    auto mir_2 = mirna::MirExp::init_from_alignment(aln2);
    CHECK(mir_2.value == 0.5);
    CHECK(mir_2.tails[2].lens[8].value == 0.5); // G tail, len 8 (10 - 2 tail len)
    CHECK(mir_2.tails[2].lens[10].value == 0); // G tail, len 10
    CHECK(mir_2.tails[0].lens[8].value == 0); // A tail, len 8
  }
}

TEST_CASE("MirExp: operator+ and partial exp") {
  auto mir_1 = mirna::MirExp {
    total_value, // value 
    {{
      a_tail,
      c_tail,
      g_tail,
      t_tail,
      o_tail,
      gm_exp
    }}
  };
  
  auto g_tail_copy = g_tail;
  auto gm_exp_copy = gm_exp;

  g_tail_copy.lens[20] = {0.0};
  g_tail_copy.lens[18] = {0.5};
  g_tail_copy.lens[21] = {0.25};
  g_tail_copy.lens[26] = {0.25};  
  gm_exp_copy.lens[18] = {1.5};
  gm_exp_copy.lens[21] = {0.5};
  gm_exp_copy.lens[22] = {0.5};
  
  auto mir_2 = mirna::MirExp {
    total_value, // value 
    {{
      a_tail,
      c_tail,
      g_tail_copy,
      t_tail,
      o_tail,
      gm_exp_copy
    }}
  };  
  SECTION("operator+") {
    auto mir_3 = mir_1 + mir_2;
    CHECK(mir_3.value == 20.0);

    // A tail
    CHECK(mir_3.tails[0].value == 2.0);
    REQUIRE(!mir_3.tails[0].lens.empty());
    CHECK(mir_3.tails[0].lens[18].value == 1.0);
    CHECK(mir_3.tails[0].lens[19].value == 1.0);  
    // C tail
    CHECK(mir_3.tails[1].value == 2.0);
    REQUIRE(!mir_3.tails[1].lens.empty());
    CHECK(mir_3.tails[1].lens[18].value == 0.5);
    CHECK(mir_3.tails[1].lens[19].value == 0.5);
    CHECK(mir_3.tails[1].lens[20].value == 1.0);
    // G tail
    CHECK(mir_3.tails[2].value == 2.0);
    REQUIRE(!mir_3.tails[2].lens.empty());
    CHECK(mir_3.tails[2].lens[18].value == 0.5);
    CHECK(mir_3.tails[2].lens[19].value == 1.0);
    CHECK(mir_3.tails[2].lens[20].value == 0.5);  
    CHECK(mir_3.tails[2].lens[21].value == 0.25);  
    CHECK(mir_3.tails[2].lens[26].value == 0.25);  
    // T tail
    CHECK(mir_3.tails[3].value == 2.0);
    REQUIRE(!mir_3.tails[3].lens.empty());
    CHECK(mir_3.tails[3].lens[20].value == 0.5);
    CHECK(mir_3.tails[3].lens[21].value == 0.5);
    CHECK(mir_3.tails[3].lens[23].value == 1.0);
    
    // Other tail
    CHECK(mir_3.tails[4].value == 2.0);
    REQUIRE(!mir_3.tails[4].lens.empty());  
    CHECK(mir_3.tails[4].lens[22].value == 2.0);  
    // Genome mathcing
    CHECK(mir_3.tails[5].value == 10.0);
    REQUIRE(!mir_3.tails[5].lens.empty());  
    CHECK(mir_3.tails[5].lens[18].value == 4.0);
    CHECK(mir_3.tails[5].lens[21].value == 0.5);
    CHECK(mir_3.tails[5].lens[22].value == 0.5);
    CHECK(mir_3.tails[5].lens[23].value == 2.0);
    CHECK(mir_3.tails[5].lens[25].value == 3.0);  
  }
  
  SECTION("operator*") {
    auto mir_2 = mir_1 * 5.0;
    CHECK(mir_2.value == 50.0);

    // A tail
    CHECK(mir_2.tails[0].value == 5.0);
    REQUIRE(!mir_2.tails[0].lens.empty());
    CHECK(mir_2.tails[0].lens[18].value == 2.5);
    CHECK(mir_2.tails[0].lens[19].value == 2.5);

    // C tail
    CHECK(mir_2.tails[1].value == 5.0);
    REQUIRE(!mir_2.tails[1].lens.empty());
    CHECK(mir_2.tails[1].lens[18].value == 1.25);
    CHECK(mir_2.tails[1].lens[19].value == 1.25);
    CHECK(mir_2.tails[1].lens[20].value == 2.5);
    
    // G tail
    CHECK(mir_2.tails[2].value == 5.0);
    REQUIRE(!mir_2.tails[2].lens.empty());
    CHECK(mir_2.tails[2].lens[19].value == 2.5);
    CHECK(mir_2.tails[2].lens[20].value == 2.5);

    // T tail
    CHECK(mir_2.tails[3].value == 5.0);
    REQUIRE(!mir_2.tails[3].lens.empty());
    CHECK(mir_2.tails[3].lens[20].value == 1.25);
    CHECK(mir_2.tails[3].lens[21].value == 1.25);
    CHECK(mir_2.tails[3].lens[23].value == 2.5);
    
    // Other tail
    CHECK(mir_2.tails[4].value == 5.0);
    REQUIRE(!mir_2.tails[4].lens.empty());  
    CHECK(mir_2.tails[4].lens[22].value == 5.0);

    // Genome mathcing
    CHECK(mir_2.tails[5].value == 25.0);
    REQUIRE(!mir_2.tails[5].lens.empty());  
    CHECK(mir_2.tails[5].lens[18].value == 12.5);
    CHECK(mir_2.tails[5].lens[23].value == 5.0);
    CHECK(mir_2.tails[5].lens[25].value == 7.5);
  }

  SECTION("get_partial_exp") {
    CHECK(mir_1.get_partial_exp() == 5.0);
    CHECK(mir_2.get_partial_exp() == 5.0);
  }
}


TEST_CASE("MirExp: transform from tail based exp to len based exp") {
  auto mir_1 = mirna::MirExp {
    total_value, // value 
    {{
      a_tail,
      c_tail,
      g_tail,
      t_tail,
      o_tail,
      gm_exp
    }}
  };

  auto&& len_based_exp = mir_1.get_len_based_exp();
  
  auto& len_18_exp = len_based_exp[18];
  auto& len_19_exp = len_based_exp[19];
  auto& len_20_exp = len_based_exp[20];
  auto& len_21_exp = len_based_exp[21];
  auto& len_22_exp = len_based_exp[22];
  auto& len_23_exp = len_based_exp[23];
  auto& len_25_exp = len_based_exp[25];
  
  auto sum = [](auto&& span){
    return std::accumulate(span.begin(), span.end(), 0.0);
  };

  CHECK(sum(len_18_exp) == 3.25);
  CHECK(sum(len_19_exp) == 1.25);
  CHECK(sum(len_20_exp) == 1.25);
  CHECK(sum(len_21_exp) == 0.25);
  CHECK(sum(len_22_exp) == 1.0);
  CHECK(sum(len_23_exp) == 1.5);
  CHECK(sum(len_25_exp) == 1.5);
  CHECK(len_18_exp[0] == 0.5);
  CHECK(len_18_exp[1] == 0.25);
  CHECK(len_18_exp[2] == 0);
  CHECK(len_18_exp[3] == 0);
  CHECK(len_18_exp[4] == 0);
  CHECK(len_18_exp[5] == 2.5);

}

TEST_CASE("use case: expression matrix") {
  
  SECTION("use case: canonical mirna expression") {
    std::cout << "TODO: utility/expression/mirna/mirna_exp.cpp requires test: use case: canonical mirna expression\n";

    /// TODO: use case: canonical mirna expression
    // mirna::ExpressionMatrix<mirna::MirExp> exp_mat;
    // exp_mat["miR92a-1-3p"] = mirna::MirExp{
  }
  SECTION("use case: isomir expression") {
    std::cout << "TODO: utility/expression/mirna/mirna_exp.cpp requires test: use case: isomir expression\n";

    /// TODO: use case: isomir expression
  }
  SECTION("use case: isomir expression with mirna seed") {
    std::cout << "TODO: utility/expression/mirna/mirna_exp.cpp requires test: use case: isomir expression with mirna seed\n";

    /// TODO: use case: isomir expression with mirna seed
  }

}

TEST_CASE("use case: sized normalization") {
  std::cout << "TODO: utility/expression/mirna/mirna_exp.cpp requires test: use case: sized normalization\n";
    
  /// TODO: use case: sized normalization
}
