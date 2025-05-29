#include <biovoltron/algo/align/tailor/tailor.hpp>
#include <biovoltron/utility/istring.hpp>
#include <biovoltron/file_io/fasta.hpp> 
#include <biovoltron/file_io/fastq.hpp> 
#include <catch.hpp>

namespace ranges = std::ranges;
using namespace biovoltron;

SCENARIO("Tailor: single match") {
  auto ref = std::vector<FastaRecord<>>{
    {"chr1", "CGATCGATCGATGCATCGATAGGGTAGCTAGCTATTAAGAGCTCTCTATGAGATGCTAGACGTATGCATGAGTCCGTATCATATGCTAGCTGAGTCGTACGTAGGGGG"}, 
    {"chr2", "TAGGTTTTAGTGATCTATAGAGAAAGAAGATCTCTCCGCGCGTATACTCGTCGGCGTCATATCGACGTATATATGCGCATCATATCGAGTCGATATCC"}, 
    {"chr3", "CGATTAGGCCGATATAGCGGCGCGCCCTCTTAGAGGGATTCGAATTAGATATATTAGGGGGTTATGCAGCATCGCTTAGCTGCCGGCGCG"}, 
    {"chr4", "GATGCTATACGATGCATACTACGATGACTAGCATCGATCGACTAGCTATATAGCTCGAGCATCGATATATGACTAGTCGTAGGAATAGGG"}, 
    {"chr5", "GGAGTAGCGATAGTAGTATGCATGACTGCAGTCATGACGTATAAGAGCGACGTTAGCAGAGCACTAGTAGTACTATAC"}
  };

  auto read = FastqRecord<>{};
  read.name = "read";
  read.qual = "!!!!!!!!!!!!!!!!!!!!";

  auto index = Index{5};
  index.make_index(ref);

  auto rc_ref = ref;
  for (auto& record : rc_ref)
    record.seq = Codec::rev_comp(record.seq);
  auto rc_index = Index{5};
  rc_index.make_index(rc_ref);

  auto tailor = Tailor{index, rc_index};
  tailor.allow_seed_mismatch = true;

  SECTION("Unqulified read: too short (< seed_len)") {
    read.seq = "GATTGTTGC";
    assert(read.seq.size() < tailor.seed_len);
    const auto aln = tailor.search(read);
    CHECK(aln.hits.empty());
  }

  SECTION("Unqulified reads: with 'N' base.") {
    read.seq = "NNNNNNAATTGATTGATTGATTGATTGTTGC";
    assert(read.seq.size() >= tailor.seed_len);
    const auto aln = tailor.search(read);
    CHECK(aln.hits.empty());
  }

  GIVEN("One mismatch at seed region") {
    read.seq = ref.front().seq.substr(2, 25);
    read.seq[4] = Codec::comp(read.seq[4]);
    WHEN("Turn off allow_seed_mismatch") {
      tailor.allow_seed_mismatch = false;
      const auto aln = tailor.search(read);
      THEN("Empty align") {
        REQUIRE(aln.hits.empty());
      }
    }
  }

  SECTION("Exact match, no tail") {
    read.seq = ref.front().seq.substr(2, 25);
    const auto aln = tailor.search(read);

    CHECK(aln.name == read.name);
    CHECK(aln.seq == read.seq);
    CHECK(aln.qual == read.qual);
    CHECK(aln.forward);
    CHECK(aln.tail_pos == -1);
    CHECK(aln.counts == 0);

    REQUIRE(aln.hits.size() == 1);
    CHECK(aln.hits.front().mismatches.empty());
    CHECK(aln.hits.front().intv == Interval{"chr1", 2, 27});
  }

  SECTION("Exact match to reverse strand, no tail") {
    read.seq = ref.front().seq.substr(2, 25);
    read.seq = Codec::rev_comp(read.seq);
    const auto aln = tailor.search(read);

    CHECK(aln.name == read.name);
    CHECK(aln.seq == read.seq);
    CHECK(aln.qual == read.qual);
    CHECK(!aln.forward);
    CHECK(aln.tail_pos == -1);
    CHECK(aln.counts == 0);

    REQUIRE(aln.hits.size() == 1);
    CHECK(aln.hits.front().mismatches.empty());
    // std::cout << aln.hits.front().intv.to_string() << "\n";
    CHECK(aln.hits.front().intv == Interval{"chr1", 2, 27, '-'});
  }

  SECTION("Tail len 1") {
    read.seq = ref.front().seq.substr(2, 26);
    read.seq.back() = Codec::comp(read.seq.back());
    const auto aln = tailor.search(read);

    CHECK(aln.name == read.name);
    CHECK(aln.seq == read.seq);
    CHECK(aln.qual == read.qual);
    CHECK(aln.forward);
    CHECK(aln.tail_pos == 25);
    CHECK(aln.counts == 0);

    REQUIRE(aln.hits.size() == 1);
    CHECK(aln.hits.front().mismatches.empty());
    CHECK(aln.hits.front().intv == Interval{"chr1", 2, 27});
  }

  SECTION("One mismatch at non-seed region, no tail") {
    read.seq = ref.front().seq.substr(2, 25);
    read.seq[20] = Codec::comp(read.seq[20]);
    const auto aln = tailor.search(read);

    CHECK(aln.name == read.name);
    CHECK(aln.seq == read.seq);
    CHECK(aln.qual == read.qual);
    CHECK(aln.forward);
    CHECK(aln.tail_pos == -1);
    CHECK(aln.counts == 0);

    REQUIRE(aln.hits.size() == 1);
    REQUIRE(aln.hits.front().mismatches.size() == 1);
    CHECK(aln.hits.front().mismatches.front() == Mismatch{20, Codec::comp(read.seq[20])});
    CHECK(aln.hits.front().intv == Interval{"chr1", 2, 27});
  }

  SECTION("Tail with length 5") {
    read.seq = ref.front().seq.substr(2, 30);
    read.seq[25] = Codec::comp(read.seq[25]);
    read.seq[27] = Codec::comp(read.seq[27]);
    const auto aln = tailor.search(read);

    CHECK(aln.name == read.name);
    CHECK(aln.seq == read.seq);
    CHECK(aln.qual == read.qual);
    CHECK(aln.forward);
    CHECK(aln.tail_pos == 25);
    CHECK(aln.counts == 0);

    REQUIRE(aln.hits.size() == 1);
    CHECK(aln.hits.front().mismatches.empty());
    CHECK(aln.hits.front().intv == Interval{"chr1", 2, 27});
  }

  SECTION("One mismatch at seed region, perfect match at non-seed, no tail") {
    read.seq = ref.front().seq.substr(2, 25);
    read.seq[4] = Codec::comp(read.seq[4]);
    const auto aln = tailor.search(read);

    CHECK(aln.name == read.name);
    CHECK(aln.seq == read.seq);
    CHECK(aln.qual == read.qual);
    CHECK(aln.forward);
    CHECK(aln.tail_pos == -1);
    CHECK(aln.counts == 0);

    REQUIRE(aln.hits.size() == 1);
    REQUIRE(aln.hits.front().mismatches.size() == 1);
    CHECK(aln.hits.front().mismatches.front() == Mismatch{4, Codec::comp(read.seq[4])});
    CHECK(aln.hits.front().intv == Interval{"chr1", 2, 27});
  }

  SECTION("One mismatch at seed region, tail with length 1") {
    read.seq = ref.front().seq.substr(2, 26);
    read.seq[4] = Codec::comp(read.seq[4]);
    read.seq[25] = Codec::comp(read.seq[25]);
    const auto aln = tailor.search(read);

    CHECK(aln.name == read.name);
    CHECK(aln.seq == read.seq);
    CHECK(aln.qual == read.qual);
    CHECK(aln.forward);
    CHECK(aln.tail_pos == 25);
    CHECK(aln.counts == 0);

    REQUIRE(aln.hits.size() == 1);
    REQUIRE(aln.hits.front().mismatches.size() == 1);
    CHECK(aln.hits.front().mismatches.front() == Mismatch{4, Codec::comp(read.seq[4])});
    CHECK(aln.hits.front().intv == Interval{"chr1", 2, 27});
  }

  SECTION("One mismatch at seed region, one mismatch at non-seed region, no tail") {
    read.seq = ref.front().seq.substr(2, 25);
    read.seq[4] = Codec::comp(read.seq[4]);
    read.seq[20] = Codec::comp(read.seq[20]);
    const auto aln = tailor.search(read);

    CHECK(aln.name == read.name);
    CHECK(aln.seq == read.seq);
    CHECK(aln.qual == read.qual);
    CHECK(aln.forward);
    CHECK(aln.tail_pos == -1);
    CHECK(aln.counts == 0);

    REQUIRE(aln.hits.size() == 1);
    REQUIRE(aln.hits.front().mismatches.size() == 2);
    const auto& mms = aln.hits.front().mismatches;
    CHECK(ranges::find(mms, Mismatch{4, Codec::comp(read.seq[4])}) != mms.end());
    CHECK(ranges::find(mms, Mismatch{20, Codec::comp(read.seq[20])}) != mms.end());
    CHECK(aln.hits.front().intv == Interval{"chr1", 2, 27});
  }

  SECTION("One mismatch at seed region, tail with length 5") {
    read.seq = ref.front().seq.substr(2, 30);
    read.seq[4] = Codec::comp(read.seq[4]);
    //tail
    read.seq[25] = Codec::comp(read.seq[25]);
    read.seq[27] = Codec::comp(read.seq[27]);
    const auto aln = tailor.search(read);

    CHECK(aln.name == read.name);
    CHECK(aln.seq == read.seq);
    CHECK(aln.qual == read.qual);
    CHECK(aln.forward);
    CHECK(aln.tail_pos == 25);
    CHECK(aln.counts == 0);

    REQUIRE(aln.hits.size() == 1);
    REQUIRE(aln.hits.front().mismatches.size() == 1);
    CHECK(aln.hits.front().mismatches.front() == Mismatch{4, Codec::comp(read.seq[4])});
    CHECK(aln.hits.front().intv == Interval{"chr1", 2, 27});
  }

  SECTION("Drop two mismatch at seed region") {
    read.seq = ref.front().seq.substr(2, 25);
    read.seq[4] = Codec::comp(read.seq[4]);
    read.seq[8] = Codec::comp(read.seq[8]);
    const auto aln = tailor.search(read);
    CHECK(aln.hits.empty());
  }

  SECTION("Align T2C reads, converting C back to T") {
    read.seq = ref.front().seq.substr(2, 20);
    auto first_T_itr = ranges::find(read.seq, 'T');
    auto second_T_itr = ranges::find(first_T_itr + 1, read.seq.end(), 'T');
    *second_T_itr = 'C';
    auto third_T_itr = ranges::find(second_T_itr + 1, read.seq.end(), 'T');
    *third_T_itr = 'C';

    /// first search failed
    const auto aln = tailor.search(read);
    CHECK(aln.hits.empty());
    
    tailor.enable_c2t = true;
    tailor.allow_seed_mismatch = false;
    const auto tc_aln = tailor.search(read);

    CHECK(tc_aln.name == read.name);
    CHECK(tc_aln.seq == read.seq);
    CHECK(tc_aln.qual == read.qual);
    CHECK(tc_aln.forward);
    CHECK(tc_aln.tail_pos == -1);
    CHECK(tc_aln.counts == 0);

    REQUIRE(tc_aln.hits.size() == 1);
    CHECK(tc_aln.hits.front().mismatches.size() == 0);
    //                      5    9
    // origin:        ATCGA TCGA TGCATCGATAG
    // T2C (read):    ATCGA CCGA CGCATCGATAG
    // 
    // since tailor do RC first, during searching:
    //                          9    5   
    //                CTATCGATGCA TCGA TCGAT
    //                CTATCGATGCG TCGG TCGAT
    auto t2c_site_itr = tc_aln.hits.front().tc_set.begin();
    REQUIRE(*t2c_site_itr == 5);
    REQUIRE(*(++t2c_site_itr) == 9);
    CHECK(tc_aln.hits.front().intv == Interval{"chr1", 2, 22});
  }
}

SCENARIO("Tailor: multiple matching occurs") {
  auto ref = std::vector<FastaRecord<>>{
    {"chr1", "CGATCGATCGATGCATCGATAGGGTAGCTAGCTATTAAGAGCTCTCTATGAGATGCTAGACGTATGCATGAGTCCGTATCATATGCTAGCTGAGTCGTACGTAGGGGG"}, 
    {"chr2", "TAGGTTTTAGTGATCTATAGAGAAAGAAGATCTCTCCGCGCGTATACTCGTCGGCGTCATATCGACGTATATATGCGCATCATATCGAGTCGATATCC"}, 
    {"chr3", "CGATTAGGCCGATATAGCGGCGCGCCCTCTTAGAGGGATTCGAATTAGATATATTAGGGGGTTATGCAGCATCGCTTAGCTGCCGGCGCG"}, 
    {"chr4", "GATGCTATACGATGCATACTACGATGACTAGCATCGATCGACTAGCTATATAGCTCGAGCATCGATATATGACTAGTCGTAGGAATAGGG"}, 
    {"chr5", "GGAGTAGCGATAGTAGTATGCATGACTGCAGTCATGACGTATAAGAGCGACGTTAGCAGAGCACTAGTAGTACTATAC"}
  };

  auto read = FastqRecord<>{};
  read.name = "read";
  read.qual = "!!!!!!!!!!!!!!!!!!!!";

  WHEN("Found 2 perfect match (with non-seed mismatch, but different base substitution) and 1 match with tail.") {
    std::copy_n(ref[0].seq.begin() + 2, 25, ref[1].seq.begin() + 5);
    std::copy_n(ref[0].seq.begin() + 2, 25, ref[2].seq.begin() + 6);
    read.seq = ref[0].seq.substr(2, 25);

    auto mismatch_pos = 20;
    auto correct_base1 = 'A';
    auto correct_base2 = 'T';
    auto wrong_base = 'C';

    read.seq[mismatch_pos] = wrong_base;
    ref[0].seq[mismatch_pos + 2] = correct_base1;
    ref[1].seq[mismatch_pos + 5] = correct_base2;
    // with tail
    ref[2].seq[mismatch_pos + 6] = correct_base2;
    ref[2].seq[mismatch_pos + 6 + 3] = Codec::comp(read.seq[mismatch_pos + 3]);

    auto index = Index{5};
    index.make_index(ref);

    auto rc_ref = ref;
    for (auto& record : rc_ref)
      record.seq = Codec::rev_comp(record.seq);
    auto rc_index = Index{5};
    rc_index.make_index(rc_ref);

    auto tailor = Tailor{index, rc_index};
    tailor.allow_seed_mismatch = true;
    const auto aln = tailor.search(read);

    THEN("Report shorter tail") {
      CHECK(aln.name == read.name);
      CHECK(aln.seq == read.seq);
      CHECK(aln.qual == read.qual);
      CHECK(aln.forward);
      CHECK(aln.tail_pos == -1);
      CHECK(aln.counts == 0);

      CHECK(aln.hits.size() == 2);
      for (const auto& hit : aln.hits) {
        REQUIRE(hit.mismatches.size() == 1);
        CHECK((
          (hit.mismatches[0] == Mismatch{mismatch_pos, correct_base1}) || 
          (hit.mismatches[0] == Mismatch{mismatch_pos, correct_base2}) 
        ));
        CHECK((
          (hit.intv == Interval{"chr1", 2, 27}) ||
          (hit.intv == Interval{"chr2", 5, 30})
        ));
      }
    }
  }

  WHEN("Report same tail length(0) and same number of mismatch") {
    std::copy_n(ref[0].seq.begin() + 2, 25, ref[1].seq.begin() + 5);
    std::copy_n(ref[0].seq.begin() + 2, 25, ref[2].seq.begin() + 6);
    read.seq = ref[0].seq.substr(2, 25);

    auto mismatch_pos1 = 3;
    auto mismatch_pos2 = 4;
    auto mismatch_pos3 = 5;

    ref[0].seq[mismatch_pos1 + 2] = Codec::comp(read.seq[mismatch_pos1]);
    ref[1].seq[mismatch_pos2 + 5] = Codec::comp(read.seq[mismatch_pos2]);
    ref[2].seq[mismatch_pos3 + 6] = Codec::comp(read.seq[mismatch_pos3]);
    auto correct_base = ref[2].seq[mismatch_pos3 + 6];

    auto index = Index{5};
    index.make_index(ref);

    auto rc_ref = ref;
    for (auto& record : rc_ref)
      record.seq = Codec::rev_comp(record.seq);
    auto rc_index = Index{5};
    rc_index.make_index(rc_ref);

    auto tailor = Tailor{index, rc_index};
    tailor.allow_seed_mismatch = true;
    const auto aln = tailor.search(read);

    THEN("Report the one with mismatch toward 3'") {
      CHECK(aln.name == read.name);
      CHECK(aln.seq == read.seq);
      CHECK(aln.qual == read.qual);
      CHECK(aln.forward);
      CHECK(aln.tail_pos == -1);
      CHECK(aln.counts == 0);

      REQUIRE(aln.hits.size() == 1);
      REQUIRE(aln.hits[0].mismatches.size() == 1);
      CHECK(aln.hits[0].mismatches[0] == Mismatch{mismatch_pos3, correct_base});
      CHECK(aln.hits[0].intv == Interval{"chr3", 6, 31});
    }
  }

  GIVEN("Found matches on both forward and reverse") {

    WHEN("Reverse report perfect match, forward report 1 tail match") {
      read.seq = ref[0].seq.substr(0, 26);
      read.seq.back() = Codec::comp(read.seq.back());
      auto rc_read = Codec::rev_comp(read.seq);

      ranges::copy(rc_read, ref[1].seq.end() - rc_read.size());

      auto index = Index{5};
      index.make_index(ref);

      auto rc_ref = ref;
      for (auto& record : rc_ref)
        record.seq = Codec::rev_comp(record.seq);
      auto rc_index = Index{5};
      rc_index.make_index(rc_ref);

      auto tailor = Tailor{index, rc_index};
      tailor.allow_seed_mismatch = true;
      const auto aln = tailor.search(read);

      THEN("Report shorter tail: reverse") {
        CHECK(aln.name == read.name);
        CHECK(aln.seq == read.seq);
        CHECK(aln.qual == read.qual);
        CHECK_FALSE(aln.forward);
        CHECK(aln.tail_pos == -1);
        CHECK(aln.counts == 0);

        REQUIRE(aln.hits.size() == 1);
        CHECK(aln.hits[0].mismatches.empty());
        CHECK(aln.hits[0].intv == Interval{"chr2", ref[1].seq.size() - rc_read.size(), ref[1].seq.size(), '-'});
      }
    }

    WHEN("Forward and reverse report same tail length") {
      read.seq = ref[0].seq.substr(0, 26);
      read.seq[3] = Codec::comp(read.seq[3]); // mismatch
      auto rc_read = Codec::rev_comp(read.seq);
      ranges::copy(rc_read, ref[1].seq.end() - rc_read.size());
      read.seq.back() = Codec::comp(read.seq.back());

      auto index = Index{5};
      index.make_index(ref);

      auto rc_ref = ref;
      for (auto& record : rc_ref)
        record.seq = Codec::rev_comp(record.seq);
      auto rc_index = Index{5};
      rc_index.make_index(rc_ref);

      auto tailor = Tailor{index, rc_index};
      tailor.allow_seed_mismatch = true;
      const auto aln = tailor.search(read);

      THEN("Report fewer mismatch") {
        CHECK(aln.name == read.name);
        CHECK(aln.seq == read.seq);
        CHECK(aln.qual == read.qual);
        CHECK_FALSE(aln.forward);
        CHECK(aln.tail_pos == read.seq.size() - 1);
        CHECK(aln.counts == 0);

        REQUIRE(aln.hits.size() == 1);
        CHECK(aln.hits[0].mismatches.empty());
        CHECK(aln.hits[0].intv == Interval{"chr2", ref[1].seq.size() - rc_read.size() + 1, ref[1].seq.size(), '-'});
      }
    }

    WHEN("Forward and reverse report same tail length and same number of mismatch") {
      read.seq = ref[0].seq.substr(0, 26);

      read.seq[3] = Codec::comp(read.seq[3]);
      read.seq[5] = Codec::comp(read.seq[5]);
      auto rc_read = Codec::rev_comp(read.seq);
      ranges::copy(rc_read, ref[1].seq.end() - rc_read.size());
      read.seq[5] = Codec::comp(read.seq[5]);

      auto index = Index{5};
      index.make_index(ref);

      auto rc_ref = ref;
      for (auto& record : rc_ref)
        record.seq = Codec::rev_comp(record.seq);
      auto rc_index = Index{5};
      rc_index.make_index(rc_ref);

      auto tailor = Tailor{index, rc_index};
      tailor.allow_seed_mismatch = true;
      const auto aln = tailor.search(read);

      THEN("Report the one with mismatch toward 3'") {
        CHECK(aln.name == read.name);
        CHECK(aln.seq == read.seq);
        CHECK(aln.qual == read.qual);
        CHECK_FALSE(aln.forward);
        CHECK(aln.tail_pos == -1);
        CHECK(aln.counts == 0);

        REQUIRE(aln.hits.size() == 1);
        REQUIRE(aln.hits[0].mismatches.size() == 1);
        CHECK(aln.hits[0].intv == Interval{"chr2", ref[1].seq.size() - rc_read.size(), ref[1].seq.size(), '-'});
      }
    }

    WHEN("Forward and reverse report same tail length, same number of mismatch, and same position of mismatch") {
      read.seq = ref[0].seq.substr(0, 25);
      auto rc_read = Codec::rev_comp(read.seq);
      ranges::copy(rc_read, ref[1].seq.end() - rc_read.size());

      auto index = Index{5};
      index.make_index(ref);

      auto rc_ref = ref;
      for (auto& record : rc_ref)
        record.seq = Codec::rev_comp(record.seq);
      auto rc_index = Index{5};
      rc_index.make_index(rc_ref);

      auto tailor = Tailor{index, rc_index};
      tailor.allow_seed_mismatch = true;
      const auto aln = tailor.search(read);

      THEN("Choose forward over reverse") {
        CHECK(aln.name == read.name);
        CHECK(aln.seq == read.seq);
        CHECK(aln.qual == read.qual);
        CHECK(aln.forward);
        CHECK(aln.tail_pos == -1);
        CHECK(aln.counts == 0);

        REQUIRE(aln.hits.size() == 1);
        CHECK(aln.hits[0].mismatches.empty());
        CHECK(aln.hits[0].intv == Interval{"chr1", 0, read.seq.size()});
      }
    }
  }

  WHEN("Found too many matches") {
    for (auto i = 1; i < ref.size(); i++) {
      std::copy_n(ref[0].seq.begin(), 25, ref[i].seq.begin());
      std::copy_n(ref[0].seq.begin(), 25, ref[i].seq.end() - 25);
    }
    read.seq = ref[0].seq.substr(0, 25);

    auto index = Index{5};
    index.make_index(ref);

    auto rc_ref = ref;
    for (auto& record : rc_ref)
      record.seq = Codec::rev_comp(record.seq);
    auto rc_index = Index{5};
    rc_index.make_index(rc_ref);

    auto tailor = Tailor{index, rc_index};
    tailor.allow_seed_mismatch = true;
    tailor.max_multi = 5;
    const auto aln = tailor.search(read);

    THEN("Drop this read") {
      CHECK(aln.hits.empty());
    }
  }
}
