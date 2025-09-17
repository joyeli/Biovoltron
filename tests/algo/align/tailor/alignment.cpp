#include <biovoltron/algo/align/tailor/alignment.hpp>
#include <biovoltron/utility/istring.hpp>
#include <fstream>
#include <catch.hpp>

using namespace biovoltron;
namespace ranges = std::ranges;

TEST_CASE("Alignment::aln_to_sam_list - Converts Alignment to SAM records", "[Alignment]") {

  SECTION("A list of aligments to SAM file") {
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
    auto aln = Alignment{
      "seq1",
      "AACCGGTTGG",
      "!!!!!!!!!!",
      true, // forward?
      8, // tail position
      {hit1, hit2},
      6
    };

    auto alns = std::vector<Alignment>(10, aln); 
    auto ofs = std::ofstream{"gg"};
    for (const auto& aln : alns)
      for (const auto& record : aln_to_sam_list(aln))
        ofs << record << '\n';
  }

  SECTION("Default constructor") {
    auto aln = Alignment{};
    CHECK(aln.hits.empty());
    CHECK(aln.counts == 0);
  }

  SECTION("Unmappable alignment output empyt sam list") {
    auto aln = Alignment{
      "seq1",
      "AACCGGTTGG",
      "!!!!!!!!!!",
      true, // forward?
      -1, // tail position
      {},
      6
    };

    auto sams = aln_to_sam_list(aln);
    REQUIRE(sams.empty());
  }

  SECTION("Mapping on forward strand") {

    SECTION("Unique mapping, no tail, no mismatch") {
      auto hit = Hit{{}, {}, {"chr1", 0, 10}};
      auto aln = Alignment{
        "seq1",
        "AACCGGTTGG",
        "!!!!!!!!!!",
        true, // forward?
        -1, // tail position
        {hit},
        6
      };

      auto sams = aln_to_sam_list(aln);

      REQUIRE(sams.size() == 1);

      CHECK(sams[0].qname == aln.name);
      CHECK(sams[0].flag == 0);
      CHECK(sams[0].rname == "chr1");
      CHECK(sams[0].pos == 1);
      CHECK(sams[0].mapq == 255);
      CHECK(sams[0].cigar == "10M");
      CHECK(sams[0].rnext == "*");
      CHECK(sams[0].pnext == 0);
      CHECK(sams[0].tlen == 0);
      CHECK(sams[0].seq == aln.seq);
      CHECK(sams[0].qual == aln.qual);

      REQUIRE(sams[0].optionals.size() == 1);
      CHECK(sams[0].optionals.front() == "NH:i:1");
    }

    SECTION("Multi mapping with tail and mismatch") {
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
      auto aln = Alignment{
        "seq1",
        "AACCGGTTGG",
        "!!!!!!!!!!",
        true, // forward?
        8, // tail position
        {hit1, hit2},
        6
      };
      auto tail_len = (aln.tail_pos == -1) ? 0 
        : aln.seq.size() - aln.tail_pos;

      auto sams = aln_to_sam_list(aln);

      REQUIRE(sams.size() == aln.hits.size());
      for (const auto& sam : sams) {
        CHECK(sam.qname == aln.name);
        CHECK(sam.flag == 0);
        CHECK(sam.mapq == 255 - tail_len);
        CHECK(sam.cigar == "8M2S");
        CHECK(sam.rnext == "*");
        CHECK(sam.pnext == 0);
        CHECK(sam.tlen == 0);
        CHECK(sam.seq == aln.seq);
        CHECK(sam.qual == aln.qual);

        CHECK(ranges::find(sam.optionals, "NH:i:2") != sam.optionals.end());
        CHECK(ranges::find(sam.optionals, "TL:Z:GG") != sam.optionals.end());
        CHECK(ranges::find(sam.optionals, "MD:Z:1C2T3") != sam.optionals.end());
      }

      CHECK(sams[0].rname == "chr1");
      CHECK(sams[0].pos == 1);

      CHECK(sams[1].rname == "chr2");
      CHECK(sams[1].pos == 11);
    }

    SECTION("Multi mapping, no tail, with different mismatch substitution") {
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
      auto aln = Alignment{
        "seq1",
        "AACCGGTTGG",
        "!!!!!!!!!!",
        true, // forward?
        -1, // tail position
        {hit1, hit2, hit3},
        6
      };
      auto tail_len = (aln.tail_pos == -1) ? 0 
        : aln.seq.size() - aln.tail_pos;

      auto sams = aln_to_sam_list(aln);

      REQUIRE(sams.size() == aln.hits.size());
      for (const auto& sam : sams) {
        CHECK(sam.qname == aln.name);
        CHECK(sam.flag == 0);
        CHECK(sam.mapq == 255 - tail_len);
        CHECK(sam.cigar == "10M");
        CHECK(sam.rnext == "*");
        CHECK(sam.pnext == 0);
        CHECK(sam.tlen == 0);
        CHECK(sam.seq == aln.seq);
        CHECK(sam.qual == aln.qual);

        REQUIRE(sam.optionals.size() == 2);
        CHECK(ranges::find(sam.optionals, "NH:i:3") != sam.optionals.end());
      }

      CHECK(sams[0].rname == "chr1");
      CHECK(sams[0].pos == 1);
      CHECK(ranges::find(sams[0].optionals, "MD:Z:1C2T5") != sams[0].optionals.end());

      CHECK(sams[1].rname == "chr2");
      CHECK(sams[1].pos == 11);
      CHECK(ranges::find(sams[1].optionals, "MD:Z:1C2T5") != sams[1].optionals.end());

      CHECK(sams[2].rname == "chr3");
      CHECK(sams[2].pos == 21);
      CHECK(ranges::find(sams[2].optionals, "MD:Z:1G2A5") != sams[2].optionals.end());
    }
  }

  SECTION("Mapping on reverse strand") {

    SECTION("Unique mapping, with tail and mismatch") {
      auto hit = Hit{
        {{1, 'C'}},
        {},
        {"chr1", 0, 10}
      };
      auto aln = Alignment{
        "seq1",
        "AACCGGTTGA",
        "abcdefghij",
        false, // forward?
        8, // tail position
        {hit},
        6
      };
      auto tail_len = (aln.tail_pos == -1) ? 0 
        : aln.seq.size() - aln.tail_pos;

      auto sams = aln_to_sam_list(aln);

      REQUIRE(sams.size() == 1);

      CHECK(sams[0].qname == aln.name);
      CHECK(sams[0].flag == 16);
      CHECK(sams[0].rname == "chr1");
      CHECK(sams[0].pos == 1);
      CHECK(sams[0].mapq == 255 - tail_len);
      CHECK(sams[0].cigar == "2S8M");
      CHECK(sams[0].rnext == "*");
      CHECK(sams[0].pnext == 0);
      CHECK(sams[0].tlen == 0);
      CHECK(sams[0].seq == Codec::rev_comp(aln.seq));
      CHECK(sams[0].qual == "jihgfedcba");

      CHECK(ranges::find(sams[0].optionals, "NH:i:1") != sams[0].optionals.end());
      CHECK(ranges::find(sams[0].optionals, "TL:Z:TC") != sams[0].optionals.end());
      CHECK(ranges::find(sams[0].optionals, "MD:Z:6C1") != sams[0].optionals.end());
    }
  }

  SECTION("Mapping with T2C reads") {

    auto hit = Hit{
      {{1, 'C'}},
      {2, 3},
      {"chr1", 0, 10}
    };
    auto aln = Alignment{
      "seq1",
      "AACCGGTTGA",
      "abcdefghij",
      false, // forward?
      8, // tail position
      {hit},
      6
    };
    auto tail_len = (aln.tail_pos == -1) ? 0 
      : aln.seq.size() - aln.tail_pos;

    auto sams = aln_to_sam_list(aln);

    REQUIRE(sams.size() == 1);

    CHECK(sams[0].qname == aln.name);
    CHECK(sams[0].flag == 16);
    CHECK(sams[0].rname == "chr1");
    CHECK(sams[0].pos == 1);
    CHECK(sams[0].mapq == 255 - tail_len);
    CHECK(sams[0].cigar == "2S8M");
    CHECK(sams[0].rnext == "*");
    CHECK(sams[0].pnext == 0);
    CHECK(sams[0].tlen == 0);
    CHECK(sams[0].seq == Codec::rev_comp(aln.seq));
    CHECK(sams[0].qual == "jihgfedcba");

    CHECK(ranges::find(sams[0].optionals, "NH:i:1") != sams[0].optionals.end());
    CHECK(ranges::find(sams[0].optionals, "TL:Z:TC") != sams[0].optionals.end());
    CHECK(ranges::find(sams[0].optionals, "MD:Z:6C1") != sams[0].optionals.end());
  
  }
}
