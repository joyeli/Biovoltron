#include <catch.hpp>
#include <biovoltron/algo/align/inexact_match/simd_alignment.hpp>
#include <biovoltron/utility/istring.hpp>
#include <experimental/random>

using namespace biovoltron;

TEST_CASE("SIMDAlignment::Alignment - Performs local and global alignments", "[SIMDAlignment]") {
    auto gen_dna_seq = [](int len) {
        auto seq = std::string{};
        while (len--)
            seq += "ATGC"[std::experimental::randint(0, 3)];
    return Codec::to_istring(seq);
    };

    SECTION("Local Alignment: all match") {
        int t = 1000; while (t--) {
            auto ref = gen_dna_seq(200);
            auto que = ref;

            auto aligner = SIMDAlignment{};
            auto [score, ref_begin, ref_end, que_begin, que_end, cigar] = aligner.local_align(ref, que);

            REQUIRE(score == ref.size());
            REQUIRE(ref_begin == 0);
            REQUIRE(ref_end == ref.size());
            REQUIRE(que_begin == 0);
            REQUIRE(que_end == que.size());
            REQUIRE(cigar == Cigar(std::to_string(ref.size()) + "M"));
        }
    }

    SECTION("Local Alignment: all match at some position of reference") {
        int t = 1000; while (t--) {
            auto que = gen_dna_seq(150);
            auto ref = gen_dna_seq(100) + que + gen_dna_seq(100);

            auto aligner = SIMDAlignment{};
            auto [score, ref_begin, ref_end, que_begin, que_end, cigar] = aligner.local_align(ref, que);

            REQUIRE(score == que.size());
            REQUIRE(que_begin == 0);
            REQUIRE(que_end == que.size());
            REQUIRE(ref_end - ref_begin == que.size());
            REQUIRE(cigar == Cigar(std::to_string(que.size()) + "M"));
            for (auto i = ref_begin; i < ref_end; i++) {
                REQUIRE(ref[i] == que[i - ref_begin]);
            }
        }
    }

    SECTION("Local Alignment: all match with softclip") {
        int t = 1000; while (t--) {
            auto ref = gen_dna_seq(100);
            auto que = gen_dna_seq(10) + ref + gen_dna_seq(10);

            auto aligner = SIMDAlignment{};
            auto [score, ref_begin, ref_end, que_begin, que_end, cigar] = aligner.local_align(ref, que);
            
            REQUIRE(score == ref.size());
            REQUIRE(ref_begin == 0);
            REQUIRE(ref_end == ref.size());
            REQUIRE(que_end - que_begin == ref.size());

            auto cigar_str = std::string{""};
            if (que_begin != 0)
                cigar_str += std::to_string(que_begin) + "S";
            cigar_str += std::to_string(ref.size()) + "M";
            if (que_end != que.size())
                cigar_str += std::to_string(que.size() - que_end) + "S";

            REQUIRE(cigar == Cigar(cigar_str));
        }
    }

    SECTION("Global Alignment: all match") {
        int t = 1000; while (t--) {
            auto ref = gen_dna_seq(200);
            auto que = ref;

            auto aligner = SIMDAlignment{};
            auto [score, cigar] = aligner.global_align(ref, que);

            REQUIRE(score == ref.size());
            REQUIRE(cigar == Cigar(std::to_string(ref.size()) + "M"));
        }
    }
}