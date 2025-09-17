#include <biovoltron/algo/align/mapq/mapq.hpp>
#include <catch.hpp>

using namespace biovoltron;

// Mapping quality score (Phred-scaled).
// The one you see in SAM file, -10log(10)Pr, where Pr is the
// probability that the mapping position is wrong, rounded to
// the nearest integer. A value 255 indicates that mapping
// quality is unavailable.
TEST_CASE("MAPQ::Calculation - Calculates MAPQ scores", "[MAPQ]") {
  SECTION("Basic use (integration test)") {
    auto aln1_scores = std::vector<int>{10, 9, 8, 7, 6, 5, 4, 3, 2, 1};
    auto [opt_score1, sub_score1, sub_cnt1] = get_opt_subopt_count(aln1_scores);

    auto aln2_scores = std::vector<int>{20, 9, 8, 7, 6, 5, 4, 3, 2, 2};
    auto [opt_score2, sub_score2, sub_cnt2] = get_opt_subopt_count(aln2_scores);

    auto aln_pair_scores = std::vector<int>{10, 9, 8, 7, 6, 5, 4, 3, 2, 1};
    auto [opt_score, sub_score, sub_cnt]
      = get_opt_subopt_count(aln_pair_scores);

    auto UNPAIRED_PENALTY = 2;
    auto score_unpaired = opt_score1 + opt_score2 - UNPAIRED_PENALTY;
    auto paring_success = (opt_score > score_unpaired);

    // struct MemAln {
    //   int score{};      // local best SW score
    //   int score2{};     // local 2nd SW score
    //   int sub_score{};  // global 2nd SW score
    //   int align_len{};  // alignment length
    //   int sub_n{};      // number of global 2nd SW score hits
    //   float frac_rep{};
    // };
    auto aln1 = MemAln{10, 8, 9, 20, 0, 0.5};
    auto aln2 = MemAln{10, 8, 9, 20, 0, 0.5};

    if (!paring_success) {
      auto mapq1 = mem_approx_mapq_se(aln1);
      CHECK((mapq1 >= 0 && mapq1 <= 60));

      auto mapq2 = mem_approx_mapq_se(aln2);
      CHECK((mapq2 >= 0 && mapq2 <= 60));
    } else {
      auto [mapq1, mapq2] = mem_mapq_pe(aln1, aln2, score_unpaired, opt_score,
                                        sub_score, sub_cnt);
      CHECK((mapq1 >= 0 && mapq1 <= 60));
      CHECK((mapq2 >= 0 && mapq2 <= 60));
    }
  }

  SECTION("Get optimal score, suboptimal score, suboptimal score count") {
    auto scores = std::vector<int>{10, 9, 8, 7, 6, 5, 4, 3, 2, 1};
    auto diff = 6;
    auto [opt_score, sub_score, sub_cnt] = get_opt_subopt_count(scores, diff);

    CHECK(opt_score == 10);
    CHECK(sub_score == 9);
    // The number score which >= min_score excluding optimal score
    // (min score = sub_score - diff).
    CHECK(sub_cnt == 7);
  }

  SECTION("MEM approximation mapping score for single-end read") {
    SECTION("sub_n == 0") {
      auto aln = MemAln{10, 8, 9, 20, 0, 0.5};
      auto mapq = mem_approx_mapq_se(aln);
      CHECK((mapq >= 0 && mapq <= 60));
    }

    SECTION("sub_n > 0") {
      auto aln = MemAln{10, 8, 9, 20, 5, 0.5};
      auto mapq = mem_approx_mapq_se(aln);
      CHECK((mapq >= 0 && mapq <= 60));
    }

    SECTION("Mapq is 0 when max(score2, sub_score) > score") {
      SECTION("sub_score is 0") {
        auto aln = MemAln{10, 8, 0, 20, 5, 0.5};
        auto mapq = mem_approx_mapq_se(aln);
        CHECK(mapq == 0);
      }

      SECTION("sub_score is not 0") {
        auto aln = MemAln{10, 8, 14, 20, 5, 0.5};
        auto mapq = mem_approx_mapq_se(aln);
        CHECK(mapq == 0);
      }
    }
  }

  SECTION("MEM mapping score for paired-end read") {
    auto opt_score = 10;
    auto sub_score = 9;
    auto sub_cnt = 7;
    auto score_un = 8;  // Unpaired score.

    auto aln1 = MemAln{10, 8, 9, 20, 0, 0.5};
    auto aln2 = MemAln{10, 8, 9, 20, 0, 0.5};
    auto [mapq1, mapq2]
      = mem_mapq_pe(aln1, aln2, score_un, opt_score, sub_score, sub_cnt);
    CHECK((mapq1 >= 0 && mapq1 <= 60));
    CHECK((mapq2 >= 0 && mapq2 <= 60));
  }
}
