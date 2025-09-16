#ifdef NDEBUG
#undef NDEBUG
#endif
#ifdef __AVX2__
#undef __AVX2__
#endif
#include <biovoltron/algo/align/spoa/simd_alignment_engine.hpp>
#include <iostream> //debug
#include <catch.hpp>

TEST_CASE("SIMD POA") {
    using namespace biovoltron;
    SECTION("Alignment") {
        {
            std::vector<std::string> sequences = {
                "CTTTTC", "CTATATATC"
            };
            auto alignment_engine = SimdAlignmentEngine::Create(
                    AlignmentType::kNW, 5, -6, -10, -2);
            Graph graph{};

            for (const auto& it : sequences) {
                auto alignment = alignment_engine->Align(it, graph);
                graph.AddAlignment(alignment, it);
            }
            auto msa = graph.GenerateMultipleSequenceAlignment();
            // REQUIRE(msa[0] == "CT---TTTC"); // due to issue #98
            REQUIRE(msa[1] == "CTATATATC");
        }
        {
            std::vector<std::string> sequences = {
                "CTTTTC", "CTATATATC"
            };
            auto alignment_engine = SimdAlignmentEngine::Create(
                    AlignmentType::kNW, 5, -6, -10);
            Graph graph{};

            for (const auto& it : sequences) {
                auto alignment = alignment_engine->Align(it, graph);
                graph.AddAlignment(alignment, it);
            }
            auto msa = graph.GenerateMultipleSequenceAlignment();
            REQUIRE(msa[0] == "CT-T-T-TC");
            REQUIRE(msa[1] == "CTATATATC");
        }
        {
            std::vector<std::string> sequences = {
                "ATCCTGG","ATCGCTG" 
            };
            auto alignment_engine = SimdAlignmentEngine::Create(
                    AlignmentType::kNW, 5, -1, -10);
            Graph graph{};

            for (const auto& it : sequences) {
                auto alignment = alignment_engine->Align(it, graph);
                graph.AddAlignment(alignment, it);
            }
            auto msa = graph.GenerateMultipleSequenceAlignment();
            REQUIRE(msa[0] == "ATCCTGG");
            REQUIRE(msa[1] == "ATCGCTG");
        }
        {
            std::vector<std::string> sequences = {
                "ATCCTGG","ATCGCTG" 
            };
            auto alignment_engine = SimdAlignmentEngine::Create(
                    AlignmentType::kNW, 5, -10, -1);
            Graph graph{};

            for (const auto& it : sequences) {
                auto alignment = alignment_engine->Align(it, graph);
                graph.AddAlignment(alignment, it);
            }
            auto msa = graph.GenerateMultipleSequenceAlignment();
            REQUIRE(msa[0] == "ATC-CTGG");
            REQUIRE(msa[1] == "ATCGCT-G");
        }
    }
}