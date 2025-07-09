#include <biovoltron/applications/burrow_wheeler_aligner/burrow_wheeler_aligner.hpp>
#include <iostream> //debug
#include <catch.hpp>

using namespace biovoltron;

// BWA-MEM: Seeding alignment with maximal exact match (MEM),
// and then extending seed with the affine-gap Smith-Waterman
// algorithm (SW).
//
// Note: Be aware, this aligner is specifically for
// hs37d5 dataset. It is not recommended to use it for
// other dataset.

TEST_CASE("Burrow Wheeler Aligner") {
    auto aligner = BurrowWheelerAligner{};

    std::cout << "TODO: applications/burrow_wheeler_aligner/burrow_wheeler_aligner.hpp requires test" << std::endl;
    
    /*
        SECTION("Compute single-end map quality score") {
        auto alns = std::vector<BurrowWheelerAligner::Aln>{
            {90, 80, 70, true, 60, 50, 40, 30, 20, 10, false, "30M"},
            {90, 70, 70, true, 60, 50, 40, 30, 20, 10, false, "30M"},
            {90, 60, 70, true, 60, 50, 40, 30, 20, 10, false, "30M"},
            {90, 50, 70, true, 60, 50, 40, 30, 20, 10, false, "30M"},
            {90, 40, 70, true, 60, 50, 40, 30, 20, 10, false, "30M"},
            {90, 30, 70, true, 60, 50, 40, 30, 20, 10, false, "30M"}};
        auto frac_rep = 0.5;
        auto [mapq, sub_score] = BurrowWheelerAligner::compute_se_mapq(alns, frac_rep);
        CHECK((mapq >= 0 && mapq <= 60));
        CHECK(sub_score == 70);
    }

    SECTION("Set CIGAR string of alignment result for read") {
        const auto read = Codec::to_istring("AAAGGTTAAGGTTAAGGTTAAGGTTAAGGTAAAAA");
        auto aln = BurrowWheelerAligner::Aln{};
        auto profile = s_profile{};

        SECTION("Alignment CIGAR not empty: no ops") {
            auto cigar = "2M2I2D2M";
            aln.cigar = cigar;
            aligner.set_cigar(aln, read, profile);
            CHECK(aln.cigar == cigar);
        }

        SECTION("Perfect match When the align score is same as the read size") {
            aln.score = read.size();
            aligner.set_cigar(aln, read, profile);
            CHECK(aln.cigar == std::to_string(read.size()) + 'M');
            CHECK(aln.align_len == read.size());
        }

        SECTION("Regular case: use smithwaterman to get cigar string") {
            // Note: hs37d5 usually start with a lot of 'N'.
            auto head = Codec::to_istring(
                "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
                "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
                "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN");
            auto tail = Codec::to_istring("TTTT");
            aligner.ref.seq = head + read + tail;
            assert(aligner.ref_.size() > aligner.EXTEND + read.size());

            assert(profile.read == nullptr);
            aln.ref_end = aligner.ref.seq.size();

            aligner.set_cigar(aln, read, profile);
            CHECK(aln.pos == head.size());
            CHECK(aln.score == read.size());
            CHECK(aln.cigar == std::to_string(read.size()) + "M");
            CHECK(aln.align_len == read.size());
        }
    }

    SECTION("Get best alignments") {
        auto alns = std::vector<BurrowWheelerAligner::Aln>{
            {90, 80, 70, true, 60, 50, 40, 30, 20, 10, false, "30M"},
            {70, 75, 70, true, 60, 50, 40, 30, 20, 10, false, "30M"},
            {70, 75, 70, true, 60, 50, 40, 30, 20, 10, false, "30M"},
            {70, 75, 70, true, 60, 50, 40, 30, 20, 10, false, "30M"},
            {50, 60, 70, true, 60, 50, 40, 30, 20, 10, false, "30M"},
            {40, 50, 70, true, 60, 50, 40, 30, 20, 10, false, "30M"},
            {30, 40, 70, true, 60, 50, 40, 30, 20, 10, false, "30M"},
            {20, 30, 70, true, 60, 50, 40, 30, 20, 10, false, "30M"}};
        // BurrowWheelerAligner::finalize_alns(alns);

        const auto read = Codec::to_istring("AGTGACAGTGACAGTGACAGTGACAGTGAC");
        const auto rread = Codec::rev_comp(read);
        auto profile = s_profile{};
        auto rprofile = s_profile{};
        const auto frac_rep = 0.5;

        const auto aln
            = aligner.get_best_one(alns, read, rread, profile, rprofile, frac_rep);

        CHECK(aln.pos == alns.front().pos);
        CHECK(aln.score == alns.front().score);
        CHECK(aln.score2 == alns.front().score2);
        CHECK(aln.forward == alns.front().forward);
        CHECK(aln.read_end == alns.front().read_end);
        CHECK(aln.ref_end == alns.front().ref_end);
        CHECK(aln.find_cnt == alns.front().find_cnt);
        CHECK(aln.align_len == alns.front().align_len);
        CHECK(aln.rescued == alns.front().rescued);

        CHECK((aln.mapq >= 0 && aln.mapq <= 60));
        CHECK(aln.sub_score == 75);
        CHECK(aln.cigar == "30M");
    }

    SECTION("Get cigar and alignment score") {
        // beg seg: ref[:5]
        // end seg: ref[-5:]
        const auto forward = true;
        const auto beg = Codec::to_istring("AAAAA");
        const auto end = Codec::to_istring("GGGGG");
        const auto mid = Codec::to_istring("TTGGAACC");
        const auto ref = beg + mid + end;
        auto comp = [](const auto c) { return c == 4 ? 4 : 3 - c; };

        auto cmp_beg = istring{};
        std::ranges::transform(beg, std::back_inserter(cmp_beg), comp);
        auto cmp_mid = istring{};
        std::ranges::transform(mid, std::back_inserter(cmp_mid), comp);

        SECTION("Mid not equal") {
            SECTION("Beg not equal, return empty") {
                const auto read = cmp_beg + cmp_mid + end;
                const auto [score, cigar] = BurrowWheelerAligner::get_score(read, ref, forward);
                CHECK(score == 0);
                CHECK(cigar == "");
            }

            SECTION("End not equal, return empty") {
                auto cmp_end = istring{};
                std::ranges::transform(end, std::back_inserter(cmp_end), comp);

                const auto read = beg + cmp_mid + cmp_end;
                const auto [score, cigar] = BurrowWheelerAligner::get_score(read, ref, forward);
                CHECK(score == 0);
                CHECK(cigar == "");
            }

            SECTION("Beg, end equal, mid more than 1 mismatch, return empty") {
                const auto read = beg + cmp_mid + end;
                const auto [score, cigar] = BurrowWheelerAligner::get_score(read, ref, forward);
                CHECK(score == 0);
                CHECK(cigar == "");
            }

            SECTION("Beg, end equal, mid 1 mismatch") {
                auto one_mm_mid = mid;
                one_mm_mid[0] = comp(one_mm_mid[0]);
                const auto read = beg + one_mm_mid + end;
                const auto full_score = read.size();
                const auto [score, cigar] = BurrowWheelerAligner::get_score(read, ref, forward);
                CHECK(score == full_score - 5);
                CHECK(cigar == std::to_string(full_score) + "M");
            }
        }

        SECTION("Beg, mid, end equal, return fullscore") {
            const auto read = ref;
            const auto full_score = read.size();
            const auto [score, cigar] = BurrowWheelerAligner::get_score(read, ref, forward);
            CHECK(score == full_score);
            CHECK(cigar == std::to_string(full_score) + "M");
        }

        SECTION("Mid equal, begin, end not equla") {
            SECTION("Score < full score - 5, return empty") {
                auto cmp_beg = istring{};
                std::ranges::transform(beg, std::back_inserter(cmp_beg), comp);
                auto cmp_end = istring{};
                std::ranges::transform(end, std::back_inserter(cmp_end), comp);

                const auto read = cmp_beg + mid + cmp_end;
                const auto [score, cigar] = BurrowWheelerAligner::get_score(read, ref, forward);
                CHECK(score == 0);
                CHECK(cigar == "");
            }

            SECTION("Other wise, properly calculate the score") {
                auto new_beg = beg;
                new_beg[0] = Codec::comp(new_beg[0]);
                const auto read = new_beg + mid + end;
                const auto [score, cigar] = BurrowWheelerAligner::get_score(read, ref, forward);
                CHECK(score == (beg.size() - 1) + mid.size() + end.size());
                CHECK(cigar
                            == "1S"
                                     + std::to_string((beg.size() - 1) + mid.size() + end.size())
                                     + "M");
            }
        }
    }

    SECTION("paring2") {
        auto alns1 = std::vector<BurrowWheelerAligner::Aln>{
            {10, 10, 70, true, 60, 50, 40, 30, 20, 10, false, "30M"},
            {20, 20, 70, true, 60, 50, 40, 30, 20, 10, false, "30M"},
            {30, 30, 70, true, 60, 50, 40, 30, 20, 10, false, "30M"},
            {90, 40, 70, true, 60, 50, 40, 30, 20, 10, false, "30M"}};
        auto alns2 = std::vector<BurrowWheelerAligner::Aln>{
            {8, 10, 70, false, 60, 50, 40, 30, 20, 10, false, "30M"},
            {10, 90, 70, false, 60, 50, 40, 30, 20, 10, false, "30M"},
            {14, 80, 70, false, 60, 50, 40, 30, 20, 10, false, "30M"},

            {17, 70, 70, false, 60, 50, 40, 30, 20, 10, false, "30M"},
            {21, 90, 70, false, 60, 50, 40, 30, 20, 10, false, "30M"},
            {24, 10, 70, false, 60, 50, 40, 30, 20, 10, false, "30M"},

            {29, 20, 70, false, 60, 50, 40, 30, 20, 10, false, "30M"},
            {33, 40, 70, false, 60, 50, 40, 30, 20, 10, false, "30M"},
            {34, 50, 70, false, 60, 50, 40, 30, 20, 10, false, "30M"}};
        const auto PAIR_DIST = 5;
        auto aln_pairs = BurrowWheelerAligner::pairing2(alns1, alns2, PAIR_DIST);
        CHECK(std::ranges::is_sorted(aln_pairs, std::ranges::greater{},
                                                                 &BurrowWheelerAligner::AlnPair::score));
        for (const auto& [aln1, aln2] : aln_pairs) {
            CHECK(aln1.forward != aln2.forward);
            CHECK(
                (aln2.pos >= aln1.pos - PAIR_DIST && aln2.pos < aln1.pos + PAIR_DIST));
        }
    }

    SECTION(
        "Split the read, use 'N' as a delimiter. Each split have to be greater "
        "than a threshold") {
        const auto read = Codec::to_istring(
            "AAAAAAAAAAAAAAAAAAANGGGGGGGGGGGGGGGGGGGNCCCNTTTTTTTTTTTTTTTTTTT");
        const auto frags = aligner.split_read(read);
        REQUIRE(frags.size() == 3);
        for (const auto& frag : frags) CHECK(frag.size() >= aligner.SEED_LEN);
        CHECK(frags[0] == 0000000000000000000_s);
        CHECK(frags[1] == 2222222222222222222_s);
        CHECK(frags[2] == 3333333333333333333_s);
    }

    SECTION("Finalize alignments: uniquify and trim smaller score") {
        SECTION("Alignments size <= 1, no change") {
            auto alns = std::vector<BurrowWheelerAligner::Aln>{
                {90, 80, 70, true, 60, 50, 40, 30, 20, 10, false, "30M"}};
            BurrowWheelerAligner::finalize_alns(alns);
            CHECK(alns.size() == 1);
        }

        SECTION("Alignments size > 1") {
            auto alns = std::vector<BurrowWheelerAligner::Aln>{
                {90, 80, 70, true, 60, 50, 40, 30, 20, 10, false, "30M"},
                {70, 70, 70, true, 60, 50, 40, 30, 20, 10, false, "30M"},
                {70, 70, 70, true, 60, 50, 40, 30, 20, 10, false, "30M"},
                {70, 70, 70, true, 60, 50, 40, 30, 20, 10, false, "30M"},
                {50, 60, 70, true, 60, 50, 40, 30, 20, 10, false, "30M"},
                {40, 50, 70, true, 60, 50, 40, 30, 20, 10, false, "30M"},
                {30, 40, 70, true, 60, 50, 40, 30, 20, 10, false, "30M"},
                {20, 30, 70, true, 60, 50, 40, 30, 20, 10, false, "30M"}};
            BurrowWheelerAligner::finalize_alns(alns);
            REQUIRE(alns.size() == 4);
            for (const auto& aln : alns)
                REQUIRE(aln.score >= 80 - aligner.MAX_SW_DIFF);
        }
    }

    // filter out alignment score < best score - MAX_SW_DIFF
    SECTION("Filter alignments") {
        auto alns = std::vector<BurrowWheelerAligner::Aln>{
            {90, 80, 70, true, 60, 50, 40, 30, 20, 10, false, "30M"},
            {90, 70, 70, true, 60, 50, 40, 30, 20, 10, false, "30M"},
            {90, 60, 70, true, 60, 50, 40, 30, 20, 10, false, "30M"},
            {90, 50, 70, true, 60, 50, 40, 30, 20, 10, false, "30M"},
            {90, 40, 70, true, 60, 50, 40, 30, 20, 10, false, "30M"},
            {90, 30, 70, true, 60, 50, 40, 30, 20, 10, false, "30M"}};

        std::ranges::sort(alns, std::ranges::greater(), &BurrowWheelerAligner::Aln::score);

        BurrowWheelerAligner::filter_alns(alns);

        REQUIRE(alns.size() == 4);
        for (const auto& aln : alns) REQUIRE(aln.score >= 80 - aligner.MAX_SW_DIFF);
    }

    */
}




