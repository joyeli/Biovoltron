#ifdef NDEBUG
#undef NDEBUG
#endif

#include <biovoltron/algo/align/wfpoa/alignment_engine.hpp>
#include <iostream> //debug
#include <catch.hpp>

TEST_CASE("WaveFront POA") {
    using namespace biovoltron;
    SECTION("Alignment") {
        {
            auto alignment_engine = biovoltron::SimdAlignmentEngine::Create(
            biovoltron::AlignmentType::kNW, 0, -1, -1);
            biovoltron::Graph graph{};
            std::vector<std::string> sequences = {"ACTG", "ACTCG", "ATCGG"};

            for (const auto& it : sequences) {
                auto alignment = alignment_engine->Align(it, graph, 0);
                graph.AddAlignment(alignment, it);
            }

            auto c = graph.GenerateMultipleSequenceAlignment(true);
            REQUIRE(c[0] == "ACT--G");
            REQUIRE(c[1] == "ACTC-G");
            REQUIRE(c[2] == "A-TCGG");
            REQUIRE(c[3] == "ACTCGG");
        }

        {
            auto alignment_engine = biovoltron::SimdAlignmentEngine::Create(
            biovoltron::AlignmentType::kNW, 0, -1, -1);
            biovoltron::Graph graph{};
            std::vector<std::string> sequences = {"AGCTAGTGTCAATGGCTACTTTTCAGGTCCT", 
                                                "AACTAAGTGTCGGTGGCTACTATATATCAGGTCCT"};

            for (const auto& it : sequences) {
                auto alignment = alignment_engine->Align(it, graph, 0);
                graph.AddAlignment(alignment, it);
            }

            auto c = graph.GenerateMultipleSequenceAlignment(true);
            REQUIRE(c[0] == "AGCT-AGTGTCAATGGCTACT-T-T-TCAGGTCCT");
            REQUIRE(c[1] == "AACTAAGTGTCGGTGGCTACTATATATCAGGTCCT");
            REQUIRE(c[2] == "AACTAAGTGTCGGTGGCTACTATATATCAGGTCCT");
        }

        {
            auto alignment_engine = biovoltron::SimdAlignmentEngine::Create(
            biovoltron::AlignmentType::kNW, 0, -1, -1);
            biovoltron::Graph graph{};
            std::string sequence = "AG";
            std::vector<std::uint32_t> weight = {1, 2, 3};
            auto alignment = alignment_engine->Align(sequence, graph, 0);
            REQUIRE_THROWS_AS(graph.AddAlignment(alignment, sequence, weight), std::invalid_argument);
            REQUIRE_THROWS_AS(graph.AddAlignment({{1, 0}, {-1, 1}}, sequence, weight), std::invalid_argument);
            REQUIRE_THROWS_AS(graph.AddAlignment({{1, 0}, {-1, 2}}, sequence, weight), std::invalid_argument);
            REQUIRE_THROWS_AS(graph.AddAlignment({{-1, 0}, {-1, 1}}, sequence, weight), std::invalid_argument);
            REQUIRE_THROWS_AS(graph.GenerateConsensus(nullptr), std::invalid_argument);
        }
    }
}