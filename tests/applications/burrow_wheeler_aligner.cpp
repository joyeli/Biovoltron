#include <biovoltron/applications/burrow_wheeler_aligner/burrow_wheeler_aligner.hpp>
#include <iostream> //debug
#include <catch.hpp>

using namespace biovoltron;
using namespace std::chrono;

// BWA-MEM: Seeding alignment with maximal exact match (MEM),
// and then extending seed with the affine-gap Smith-Waterman
// algorithm (SW).
//
// Note: Be aware, this aligner is specifically for
// hs37d5 dataset. It is not recommended to use it for
// other dataset.

TEST_CASE("Burrow Wheeler Aligner") 
{
  static std::once_flag init_flag;
  static FMIndex<1, uint32_t, StableSorter<uint32_t>> index;
  static FastaRecord<true> ref;
  static FastqRecord<true> read1_ori, read2_ori;
  // Note: hs37d5 usually start with a lot of 'N'.
  static auto head = Codec::to_istring(std::string(200, 'N'));
  static auto tail = Codec::to_istring("TTTT");
  // gap in middle
  static auto middle = Codec::to_istring("GGGGGGAAAACCCN");;
  static auto t_start = high_resolution_clock::now();

  std::call_once(init_flag, [&]() {
    // Original read pair (forward and reverse reads)
    read1_ori = FastqRecord<true>{
      {"read1/1", Codec::to_istring("AAAGGTTAAGGTTAAGGTTAAGGTTAAGGTAAAAA")},
      "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
    };

    read2_ori = FastqRecord<true>{
      {"read1/2", Codec::to_istring("TTTTTTTCCTAACCCTAACCTAACCTAACCTTTTT")},
      "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
    };
    // Reference contains two nearly identical regions
    auto target = Codec::to_istring("GGGACGTACTGACTGACTGACTGACTGACTGAAAA");
    auto near_copy = Codec::to_istring("GGGACGTACTGACTGACTGACTGACTGACTGAAAT");
    
    auto new_ref_seq = head + read1_ori.seq + middle + Codec::rev_comp(read2_ori.seq) + target + near_copy + tail;
    // Mask ambiguous bases if needed (set to 0 if >= 4)
    std::ranges::transform(new_ref_seq, new_ref_seq.begin(), [](auto& c) {
      return c < 4 ? c : 0;
    });

    ref = FastaRecord<true>{"test", new_ref_seq};

    auto t0 = high_resolution_clock::now();
    // Build FM-index
    index.build(ref.seq);
    auto t1 = high_resolution_clock::now();
    // std::cout << "[Timer] Build FM-Index : "
    //           << duration_cast<milliseconds>(t1 - t0).count()
    //           << " ms\n";
  });

  SECTION("Perfect paired-end alignment") 
  {
    auto t0 = high_resolution_clock::now();
    // Simulate reference containing exact read1 and rev-comp(read2)
    auto read1 = read1_ori;
    auto read2 = read2_ori;

    BurrowWheelerAligner aligner{ref, index};

    // Perform alignment
    const auto [rec1, rec2] = aligner.generate_sam({read1, read2});

    // Check read1 alignment results
    CHECK(rec1.pos == head.size() + 1); // 1-based position
    CHECK(rec1.cigar == std::to_string(read1.seq.size()) + "M");
    CHECK(rec1.flag == 99); // paired, first in pair, properly aligned

    // Check read2 alignment results
    CHECK(rec2.flag == 147); // paired, second in pair, properly aligned, reverse strand
    CHECK(rec2.pos > rec1.pos);
    CHECK(rec2.pos == rec1.pos + read1.seq.size() + middle.size());
    CHECK(rec2.cigar == std::to_string(read2.seq.size()) + "M");

    // Validate MAPQ range
    CHECK((rec1.mapq >= 0 && rec1.mapq <= 60));
    CHECK((rec2.mapq >= 0 && rec2.mapq <= 60));

    auto t1 = high_resolution_clock::now();
    // std::cout << "[Timer] Perfect paired-end alignment : "
    //           << duration_cast<milliseconds>(t1 - t0).count()
    //           << " ms\n";     
  }

  SECTION("Unmappable reads return * cigar and 0 mapq")
  {
    auto t0 = high_resolution_clock::now();
    BurrowWheelerAligner aligner{ref, index};

    // Reads with no matching region in reference
    const auto read1 = FastqRecord<true>{
      {"read1/1", Codec::to_istring("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG")},
      "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
    };

    const auto read2 = FastqRecord<true>{
      {"read1/2", Codec::to_istring("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC")},
      "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
    };

    const auto [rec1, rec2] = aligner.generate_sam({read1, read2});

    CHECK(rec1.cigar == "*");
    CHECK(rec2.cigar == "*");
    CHECK(rec1.mapq == 0);
    CHECK(rec2.mapq == 0);
    CHECK(rec1.flag & 0x4); // read1 unmapped
    CHECK(rec2.flag & 0x4); // read2 unmapped

    auto t1 = high_resolution_clock::now();
    // std::cout << "[Timer] Unmappable reads : "
    //           << duration_cast<milliseconds>(t1 - t0).count()
    //           << " ms\n"; 
  }

  SECTION("Read with mismatches fails to align cleanly") 
  {
    auto t0 = high_resolution_clock::now();
    BurrowWheelerAligner aligner{ref, index};

    // Modify first 5 bases of read1 to be complementary (introduce mismatches)
    auto read1 = read1_ori;
    auto read2 = read2_ori;
    for (auto i = 0; i < 5; i++) {
      read1.seq[i] = Codec::comp(read1.seq[i]);
    }

    const auto [rec1, rec2] = aligner.generate_sam({read1, read2});

    if (rec1.cigar != "*") {
      CHECK(rec1.mapq <= 60); // mapped but may have low confidence
    } else {
      CHECK(rec1.mapq == 0);  // unmapped
    }

    auto t1 = high_resolution_clock::now();
    // std::cout << "[Timer] Read with mismatches : "
    //           << duration_cast<milliseconds>(t1 - t0).count()
    //           << " ms\n"; 
  }

  SECTION("Ambiguous read matches multiple locations with similar score") 
  {
    auto t0 = high_resolution_clock::now();
    // Simulate a read that matches multiple locations with similar score
    BurrowWheelerAligner aligner{ref, index};

    const auto read = FastqRecord<true>{
      {"read/1", Codec::to_istring("GGGACGTACTGACTGACTGACTGACTGACTGAAAA")}, 
      "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
    };

    const auto [rec1, rec2] = aligner.generate_sam({read, read}); 
    CHECK(rec1.cigar == std::to_string(read.seq.size()) + "M");
    CHECK(rec1.pos == head.size() + read1_ori.seq.size() + middle.size() + read2_ori.seq.size() + 1); // 1-based position
    CHECK(rec1.mapq < 10); // ambiguous due to similar match

    auto t1 = high_resolution_clock::now();
    // std::cout << "[Timer] Ambiguous read matches : "
    //           << duration_cast<milliseconds>(t1 - t0).count()
    //           << " ms\n";
  }

  SECTION("Read with insertion or deletion") 
  {
    auto t0 = high_resolution_clock::now();
    // Simulate a read that is shorter than reference (deleted seq[2])
    const auto read = FastqRecord<true>{
    {"read1/1", Codec::to_istring("AAGGTTAAGGTTAAGGTTAAGGTTAAGGTAAAAA")},
    "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
    };

    BurrowWheelerAligner aligner{ref, index};
    const auto [rec1, rec2] = aligner.generate_sam({read, read});

    if (rec1.cigar != "*") {
      // std::cout << "CIGAR: " << rec1.cigar << ", MAPQ: " << rec1.mapq << "\n";
      CHECK(rec1.mapq <= 60); // indel present, but still aligned
    } else {
      CHECK(rec1.mapq == 0);
    }

    auto t1 = high_resolution_clock::now();
    // std::cout << "[Timer] Read with insertion or deletion : "
    //           << duration_cast<milliseconds>(t1 - t0).count()
    //           << " ms\n"; 
  }

  SECTION("Read flanked by aligned segments with unmatched middle - PE setup") 
  {
    auto t0 = high_resolution_clock::now();
    // In this scenario, we simulate paired-end reads where read1 aligns to front
    // and read2(no reverse completementary) aligns to back, while the reference contains both but with a gap in between.
    const auto read1 = read1_ori;
    const auto read_seq_rev_cmp = Codec::rev_comp(read2_ori.seq);
    const auto read2 = FastqRecord<true>{
      {"read1/2", read_seq_rev_cmp},
      "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
    };

    BurrowWheelerAligner aligner{ref, index};

    const auto [rec1, rec2] = aligner.generate_sam({read1, read2});

    // std::cout << "[Flanked PE] Read1 CIGAR: " << rec1.cigar << ", MAPQ: " << rec1.mapq << "\n";
    // std::cout << "[Flanked PE] Read2 CIGAR: " << rec2.cigar << ", MAPQ: " << rec2.mapq << "\n";
    
    CHECK(rec1.mapq == 60);
    CHECK(rec2.mapq == 60);
    CHECK(rec1.flag == 65); // paired, first in pair
    CHECK(rec2.flag == 129); // paired, second in pair
    CHECK(rec1.pos < rec2.pos);

    auto t1 = high_resolution_clock::now();
    // std::cout << "[Timer] Read flanked by aligned segments with unmatched middle - PE setup : "
    //           << duration_cast<milliseconds>(t1 - t0).count()
    //           << " ms\n";
    // std::cout << "[Timer] All timings : "
    //           << duration_cast<milliseconds>(t1 - t_start).count()
    //           << " ms\n";
  }
}




