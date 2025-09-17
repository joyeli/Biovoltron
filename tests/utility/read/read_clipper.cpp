#include <biovoltron/utility/read/read_clipper.hpp>
#include <iostream> //debug
#include <catch.hpp>

#include <biovoltron/file_io/cigar.hpp>

using namespace biovoltron;
TEST_CASE("ReadClipper::Clipping - Clips reads based on various criteria", "[ReadUtils]") {

  SECTION("hard_clip_soft_clipped_bases") {
    SamRecord<> read;

    // case 1: front_length = 5, front_op = 'S'
    read.seq = "AAAAACCCCCGGGGG";
    read.qual = "ABCDEFGHIJKLMNO";
    read.cigar = Cigar("5S10M");

    ReadClipper::hard_clip_soft_clipped_bases(read);

    REQUIRE(read.seq == "CCCCCGGGGG");
    REQUIRE(read.qual == "FGHIJKLMNO");

    // case 2: back_length = 5, back_op = 'S'
    read.seq = "TTTTTGGGGGCCCCC";
    read.qual = "ABCDEFGHIJKLMNO";
    read.cigar = Cigar("10M5S");

    ReadClipper::hard_clip_soft_clipped_bases(read);

    REQUIRE(read.seq == "TTTTTGGGGG");
    REQUIRE(read.qual == "ABCDEFGHIJ");

    // case 3: both front and back are 'S'
    read.seq = "GGGGGTTTTTAAAAA";
    read.qual = "ABCDEFGHIJKLMNO";
    read.cigar = Cigar("5S5M5S");

    ReadClipper::hard_clip_soft_clipped_bases(read);

    REQUIRE(read.seq == "TTTTT");
    REQUIRE(read.qual == "FGHIJ");

    // case 4: no 'S' at all
    read.seq = "ACGTACGTACGT";
    read.qual = "ABCDEFGHIJKL";
    read.cigar = Cigar("12M");

    ReadClipper::hard_clip_soft_clipped_bases(read);

    REQUIRE(read.seq == "ACGTACGTACGT");
    REQUIRE(read.qual == "ABCDEFGHIJKL");
  }
  
  SECTION("revert_soft_clipped_bases") {
    // case 1: Forward strand: revert front soft-clip, adjust POS and cigar
    SamRecord<> read;
    read.flag = 0; // forward
    read.seq = "ACGTACGTAC";
    read.qual = "ABCDEFGHIJ";
    read.pos = 6; // begin() = 5
    read.cigar = Cigar("5S5M");

    ReadClipper::revert_soft_clipped_bases(read);

    REQUIRE(read.seq == "ACGTACGTAC");
    REQUIRE(read.qual == "ABCDEFGHIJ");
    REQUIRE(read.cigar.front() == Cigar::Element{5, 'M'}); // front S -> M
    REQUIRE(read.pos == 1);  // 5 - 5 + 1

    // case 2: Forward strand: no pos adjust when begin < front_length
    read.flag = 0;
    read.seq = "ACGTACGT";
    read.qual = "ABCDEFGH";
    read.pos = 2; // begin() = 1
    read.cigar = Cigar("5S3M");

    ReadClipper::revert_soft_clipped_bases(read);

    REQUIRE(read.pos == 2);
    REQUIRE(read.cigar.front() == Cigar::Element{5, 'S'});
    REQUIRE(read.seq == "ACGTACGT");
    REQUIRE(read.qual == "ABCDEFGH");

  // case 3: remove back soft-clip from seq/qual
    read.flag = 0;
    read.seq = "ACGTACGTAC";
    read.qual = "ABCDEFGHIJ";
    read.pos = 1;
    read.cigar = Cigar{"5M5S"};

    ReadClipper::revert_soft_clipped_bases(read);

    REQUIRE(read.seq == "ACGTA");
    REQUIRE(read.qual == "ABCDE");

  // case 4: Reverse strand: remove front soft-clip from seq/qual
    read.flag = SamUtil::READ_REVERSE_STRAND;
    read.seq = "TTTTTGGGGG";
    read.qual = "ABCDEFGHIJ";
    read.cigar = Cigar{"5S5M"};

    ReadClipper::revert_soft_clipped_bases(read);

    REQUIRE(read.seq == "GGGGG");
    REQUIRE(read.qual == "FGHIJ");

    // case 5: Reverse strand: back soft-clip should change to M
    read.flag = SamUtil::READ_REVERSE_STRAND;
    read.seq = "GGGGGTTTTT";
    read.qual = "ABCDEFGHIJ";
    read.cigar = Cigar{"5M5S"};

    ReadClipper::revert_soft_clipped_bases(read);

    REQUIRE(read.seq == "GGGGGTTTTT");
    REQUIRE(read.qual == "ABCDEFGHIJ");
    REQUIRE(read.cigar.back() == Cigar::Element{5, 'M'});

    // case 6: No soft-clips, should remain unchanged
    read.flag = 0;
    read.seq = "ACGT";
    read.qual = "ABCD";
    read.pos = 2;
    read.cigar = Cigar{"4M"};

    auto original = read;

    ReadClipper::revert_soft_clipped_bases(read);

    REQUIRE(read.seq == original.seq);
    REQUIRE(read.qual == original.qual);
    REQUIRE(read.pos == original.pos);
    REQUIRE(read.cigar == original.cigar);
  }
  
  SECTION("hard_clip_to_interval") {
    // case 1: no clipping needed
    SamRecord<> read;
    read.rname = "chr1";
    read.seq = "ACGTACGT";
    read.qual = "ABCDEFGH";
    read.pos = 6; // begin() = 5
    read.cigar = Cigar{"8M"}; // ref_size = 8 → end() = 13

    Interval interval{"chr1", 5, 13, '+'};

    ReadClipper::hard_clip_to_interval(read, interval);

    REQUIRE(read.seq == "ACGTACGT");
    REQUIRE(read.qual == "ABCDEFGH");

    // case 2: clip at the front
    read.rname = "chr1";
    read.seq = "ACGTACGT";
    read.qual = "ABCDEFGH";
    read.pos = 4; // begin() = 3
    read.cigar = Cigar{"8M"}; // end() = 11

    Interval interval_2{"chr1", 5, 11, '+'};

    ReadClipper::hard_clip_to_interval(read, interval_2);

    // begin = 3 < 5 → clip 2 from left
    REQUIRE(read.seq == "GTACGT");
    REQUIRE(read.qual == "CDEFGH");

    // case 3: clip at the back
    read.rname = "chr1";
    read.seq = "ACGTACGT";
    read.qual = "ABCDEFGH";
    read.pos = 5; // begin = 4
    read.cigar = Cigar{"8M"}; // end = 12

    Interval interval_3{"chr1", 4, 10, '+'};

    ReadClipper::hard_clip_to_interval(read, interval_3);

    // end = 12 > 10 → clip 2 from right
    REQUIRE(read.seq == "ACGTAC");
    REQUIRE(read.qual == "ABCDEF");

    // case 4: clip both ends
    read.rname = "chr1";
    read.seq = "ACGTACGT";
    read.qual = "ABCDEFGH";
    read.pos = 3; // begin = 2
    read.cigar = Cigar{"8M"}; // end = 10

    Interval interval_4{"chr1", 4, 8, '+'};

    ReadClipper::hard_clip_to_interval(read, interval_4);

    // left clip: 4 - 2 = 2
    // right clip: 10 - 8 = 2
    REQUIRE(read.seq == "GTAC");
    REQUIRE(read.qual == "CDEF");

    // case 5: clip size exceeds seq size
    read.rname = "chr1";
    read.seq = "ACGT";
    read.qual = "ABCD";
    read.pos = 1; // begin = 0
    read.cigar = Cigar{"4M"}; // end = 4

    Interval interval_5{"chr1", 10, 20, '+'};

    ReadClipper::hard_clip_to_interval(read, interval_5);

    // clip size = 10 - 0 = 10 → exceeds seq.size() = 4 → capped
    REQUIRE(read.seq.empty());
    REQUIRE(read.qual.empty());
  }
}

