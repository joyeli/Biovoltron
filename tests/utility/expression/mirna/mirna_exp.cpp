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
using biovoltron::Alignment;
using biovoltron::Hit;
namespace mirna = biovoltron::mirna;
///                           value, lens
auto total_value = 10.0;
auto a_tail = mirna::TailExp { 1.0, { {18, {0.5}}, {19, {0.5}} } };
auto c_tail = mirna::TailExp { 1.0, { {18, {0.25}}, {19, {0.25}}, {20, {0.5}} } };
auto g_tail = mirna::TailExp { 1.0, { {19, {0.5}}, {20, {0.5}} } };
auto t_tail = mirna::TailExp { 1.0, { {20, {0.25}}, {21, {0.25}}, {23, {0.5}} } };
auto o_tail = mirna::TailExp { 1.0, { {22, {1.0}} } };
auto gm_exp = mirna::TailExp { 5.0, { {18, {2.5}}, {23, {1.0}}, {25, {1.5}} } };

TEST_CASE("MirExp::Construction - Constructs MirExp objects", "[MirExp]") {
  
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

TEST_CASE("MirExp::Operations - Performs arithmetic and partial expression operations", "[MirExp]") {
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


TEST_CASE("MirExp::Transformation - Transforms tail-based to length-based expression", "[MirExp]") {
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

namespace{
  Hit hit_92a1{{}, {}, {"chr13", 91351361, 91351383}};
  Hit hit_92a2{{}, {}, {"chrX", 134169544, 134169566, '-'}};
  Hit hit_21_5p{{}, {}, {"chr17", 59841273, 59841295}};
  // canoncial
  Alignment r1{
    "read1", "UGAGGUAGUAGGUUGUAUAGUU", 
    std::string(22,'!'), true, -1, {hit_92a1, hit_92a2}, 1
  };

  Alignment r2{
    "read2", "UGAGGUAGUAGGUUGUAUAGUUAAA",  
    std::string(25,'!'), true, 22, {hit_92a2}, 1
  };
  // 3' nucleotide addition A-tail
  Alignment r3{
    "read3", "UGAGGUAGUAGGUUGUAUAGUUA",  
    std::string(23,'!'), true, 22, {hit_92a1}, 1
  };
  // 3' nucleotide addition U-tail
  Alignment r4{
    "read4", "UGAGGUAGUAGGUUGUAUAGUUU",
    std::string(23, '!'), true, 22, {hit_92a1, hit_92a2, hit_21_5p}, 1
  };
  // 3' upstream, 少 １ nt
  Alignment r5{
    "read5", "UGAGGUAGUAGGUUGUAUAGU",
    std::string(21, '!'), true, -1, {hit_92a1}, 1
  };   
  // 3' downstream, 長 1 nt, 模板 nt, 歸 -1 (M)
  Alignment r6{
    "read6", "UGAGGUAGUAGGUUGUAUAGUUC",  
    std::string(23,'!'), true, -1, {hit_92a1, hit_92a2}, 1
  };
  // heteropolymeric tail
  Alignment r7{
  "read7", "UGAGGUAGUAGGUUGUAUAGUUAC",
  std::string(24,'!'), true, 22, {hit_92a1}, 1
  };
  // 5' trimming
  Alignment r8{
    "read8", "GAGGUAGUAGGUUGUAUAGUU",
    std::string(21,'!'), true, -1, {hit_92a1}, 1
  };
  // 5' addition
  Alignment r9{
    "read9", "CUGAGGUAGUAGGUUGUAUAGUU",
    std::string(23,'!'), true, -1, {hit_92a1, hit_92a2}, 1
  };

  auto m1 = mirna::MirExp::init_from_alignment(r1);
  auto m2 = mirna::MirExp::init_from_alignment(r2);
  auto m3 = mirna::MirExp::init_from_alignment(r3);
  auto m4 = mirna::MirExp::init_from_alignment(r4);
  auto m5 = mirna::MirExp::init_from_alignment(r5);
  auto m6 = mirna::MirExp::init_from_alignment(r6);
  auto m7 = mirna::MirExp::init_from_alignment(r7);
  auto m8 = mirna::MirExp::init_from_alignment(r8);
  auto m9 = mirna::MirExp::init_from_alignment(r9);
  

  static auto matrix_92a() -> mirna::ExpressionMatrix<mirna::MirExp>{
    mirna::ExpressionMatrix<mirna::MirExp> exp_mat;
    exp_mat["miR92a-1-3p"] += m1;
    exp_mat["miR92a-1-3p"] += m3;
    exp_mat["miR92a-1-3p"] += m4;
    exp_mat["miR92a-1-3p"] += m5;
    exp_mat["miR92a-1-3p"] += m6;
    exp_mat["miR92a-1-3p"] += m7;
    exp_mat["miR92a-1-3p"] += m8;
    exp_mat["miR92a-1-3p"] += m9;
    return exp_mat;
  }
} // namespace

TEST_CASE("MirExp::UseCase - Expression matrix", "[MirExp]") {
  SECTION("use case: canonical mirna expression") {
    mirna::ExpressionMatrix<mirna::MirExp> exp_mat;
    exp_mat["miR92a-1-3p"] += m1;
    exp_mat["miR92a-2-3p"] += m1;
    exp_mat["miR92a-1-3p"] += m2;
    exp_mat["miR92a-1-3p"] += m3;
    exp_mat["miR-21-5p"] += m2;
    exp_mat["miR-21-5p"] += m4;

    const auto& mir92a1 = exp_mat.at("miR92a-1-3p");
    auto mir92a1_copy = mir92a1;
    auto len_base_mir92a1 = mir92a1_copy.get_len_based_exp();
    CHECK(mir92a1.value                       == Approx(2.5));
    CHECK(mir92a1.tails[0].value              == Approx(2.0));
    CHECK(mir92a1.tails[0].lens.at(22).value  == Approx(2.0));
    CHECK(mir92a1.tails[3].value              == Approx(0.0));
    CHECK(mir92a1.tails[5].value              == Approx(0.5));
    CHECK(mir92a1.tails[5].lens.at(22).value  == Approx(0.5));
    CHECK(len_base_mir92a1.at(22)[0]          == Approx(2.0));
    CHECK(len_base_mir92a1.at(22)[3]          == Approx(0.0));
    CHECK(len_base_mir92a1.at(22)[5]          == Approx(0.5));
    CHECK(mir92a1_copy.get_partial_exp()      == Approx(2.0));

    const auto& mir92a2 = exp_mat.at("miR92a-2-3p");
    auto mir92a2_copy = mir92a2;
    auto len_base_mir92a2 = mir92a2_copy.get_len_based_exp();
    CHECK(mir92a2.value                      == Approx(0.5));
    CHECK(mir92a2.tails[5].value             == Approx(0.5));
    CHECK(mir92a2.tails[5].lens.at(22).value == Approx(0.5));
    CHECK(len_base_mir92a2.at(22)[5]         == Approx(0.5));
    CHECK(mir92a2_copy.get_partial_exp()     == Approx(0.0));

    const auto& mir21_5p = exp_mat.at("miR-21-5p");
    auto mir21_5p_copy = mir21_5p;
    auto len_base_mir21_5p = mir21_5p_copy.get_len_based_exp();
    CHECK(mir21_5p.value                      == Approx(4.0/3.0));
    CHECK(mir21_5p.tails[0].value             == Approx(1.0));
    CHECK(mir21_5p.tails[3].value             == Approx(1.0/3.0));
    CHECK(mir21_5p.tails[0].lens.at(22).value == Approx(1.0));
    CHECK(mir21_5p.tails[3].lens.at(22).value == Approx(1.0/3.0));
    CHECK(len_base_mir21_5p.at(22)[0]         == Approx(1.0));
    CHECK(len_base_mir21_5p.at(22)[3]         == Approx(1.0/3.0));
    CHECK(mir21_5p_copy.get_partial_exp()     == Approx(4.0/3.0));
  }
  // 3' 
  SECTION("use case: isomir expression") {
    mirna::ExpressionMatrix<mirna::MirExp> exp_mat; 
    exp_mat["miR92a-1-3p"] += m1;
    exp_mat["miR92a-1-3p"] += m3;
    exp_mat["miR92a-1-3p"] += m4;
    exp_mat["miR92a-1-3p"] += m5;
    exp_mat["miR92a-1-3p"] += m6;
    exp_mat["miR92a-1-3p"] += m7;
    const auto& mexp = exp_mat.at("miR92a-1-3p");

    REQUIRE(!mexp.tails[0].lens.empty());
    REQUIRE(mexp.tails[2].lens.empty());
    REQUIRE(!mexp.tails[3].lens.empty());
    REQUIRE(!mexp.tails[4].lens.empty());
    REQUIRE(!mexp.tails[5].lens.empty());

    CHECK(mexp.tails[5].lens.at(22).value == Approx(0.5));   // r1
    CHECK(mexp.tails[0].lens.at(22).value == Approx(1.0));   // r3
    CHECK(mexp.tails[3].lens.at(22).value == Approx(1.0/3.0)); // r4
    CHECK(mexp.tails[5].lens.at(21).value == Approx(1.0));   // r5
    CHECK(mexp.tails[5].lens.at(23).value == Approx(0.5));   // r6
    CHECK(mexp.tails[4].lens.at(22).value == Approx(1.0)); // r7

    double sum = 0.0;
    for (int t = 0; t < (int)mexp.tails.size(); ++t){
      for(auto&& [len, lexp] : mexp.tails[t].lens){
        sum += lexp.value;
      }
    }
    CHECK(sum == Approx(mexp.value).margin(1e-12));
  }
  // 5' change seed, 3' unchange seed
  SECTION("use case: isomir expression with mirna seed") {
    auto exp_mat = matrix_92a();

    // 3' 歸到同一 seed
    // 5' (trimming, addition) 分別應歸不同 seed
    auto seed_of = [](const Alignment& aln) -> std::string{
      return aln.seq.substr(1, 7);
    };
    std::map<std::string, double> seed_bucket;
    // 加到各桶
    auto add_to_bucket = [&](const Alignment& aln) {
      seed_bucket[seed_of(aln)] += 1.0 / static_cast<double>(aln.hits.size());
    };

    add_to_bucket(r1);
    add_to_bucket(r3);
    add_to_bucket(r4);
    add_to_bucket(r5);
    add_to_bucket(r6);
    add_to_bucket(r7);
    add_to_bucket(r8);
    add_to_bucket(r9);

    const std::string seed_canonical = "GAGGUAG";   // r1
    const std::string seed_5p_trimming = "AGGUAGU"; // r8
    const std::string seed_5p_addition = "UGAGGUA"; // r9

    CHECK(seed_bucket[seed_canonical] == Approx(13.0/3.0));
    CHECK(seed_bucket[seed_5p_trimming] == Approx(1.0));
    CHECK(seed_bucket[seed_5p_addition] == Approx(0.5));

    double total = 0.0;
    for(auto& [k, v] : seed_bucket) total += v;
    const auto& mexp = exp_mat.at("miR92a-1-3p");
    CHECK(total == Approx(mexp.value));
  }

}

TEST_CASE("MirExp::UseCase - Sized normalization", "[MirExp]") {
  auto exp_mat = matrix_92a();

  // T_l = Σ(miRNA) Σ(tail) x_{l, t}
  auto len_totals = [](const mirna::ExpressionMatrix<mirna::MirExp>& mat){
    std::map<int, double> totals;
    for(auto&& [id, m] : mat){
      for(const auto& tail : m.tails){
        for(const auto& [len, lexp] : tail.lens){
          totals[len] += lexp.value;
        }
      }
    }
    return totals;
  };

  const auto totals = len_totals(exp_mat);

  // 產生每個長度的 weight
  double mean = 0.0;
  for(auto& [_, v] : totals) mean += v;
  mean /= totals.size();
  std::map<int,double> w;
  for (auto& [len, v] : totals) w[len] = (v == 0.0) ? 1.0 : (mean / v);

  // 將 w_l 乘回所有 mirna 在該長度的數值
  auto exp_norm = exp_mat;
  for (auto& [id, m] : exp_norm) {
    double total = 0.0;
    for (auto& tail : m.tails) {
      double tail_sum = 0.0;
      for (auto& [len, lexp] : tail.lens) {
        const double factor = w.count(len) ? w[len] : 1.0;
        lexp.value *= factor;
        tail_sum += lexp.value;
      }
      tail.value = tail_sum;
      total += tail_sum;
    }
    m.value = total;
  }

  // 驗等化
  const auto totals_after = len_totals(exp_norm);
  REQUIRE(totals_after.size() == totals.size());
  for (const auto& [len, val] : totals_after) {
    CHECK(val == Approx(mean).margin(1e-12));
  }
  // 驗同長度 tail 比例不變
  auto frac_at = [](const mirna::MirExp& m, int L){
    std::array<double,6> f{}; double s = 0.0;
    for (int t = 0; t < 6; ++t){
      if (auto it = m.tails[t].lens.find(L); it != m.tails[t].lens.end()){
        s += (f[t] = it->second.value);
      }
    }
    if (s > 0.0) for (auto& v : f) v /= s;
    return f;
  };
  
  const std::string mir_id = "miR92a-1-3p";
  const int L = 22;
  const auto before = frac_at(exp_mat.at(mir_id), L);
  const auto after = frac_at(exp_norm.at(mir_id), L);
  for(int t = 0; t < 6; t++) CHECK(after[t] == Approx(before[t]).margin(1e-12));
}