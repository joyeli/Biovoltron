#include <biovoltron/file_io/fastq.hpp>
#include <biovoltron/utility/istring.hpp>
#include <catch.hpp>
#include <filesystem>
#include <fstream>
#include <sstream>

using namespace biovoltron;

const auto data_path = std::filesystem::path{DATA_PATH};

inline void
CHECK_FASTQ_IDENTITY(const std::filesystem::path& path, std::string_view expected_result) {
  std::fstream fin(path, std::ios::in);
  REQUIRE(fin.is_open() == true);
  std::ostringstream oss;
  std::vector<FastqRecord<>> Records;
  FastqRecord record;
  while (fin >> record) Records.emplace_back(record);
  for (auto record : Records) oss << record << '\n';
  REQUIRE(oss.str() == expected_result);
}

TEST_CASE("FastqRecord Test", "[FastqRecord]") {
  SECTION("Default Constructor and Implicit Conversion") {
    FastqRecord<> fastqString;
    fastqString.name = "TestSequence";
    fastqString.seq = "ACGT";
    fastqString.qual = "!@?#";

    // Test implicit conversion from string to istring
    FastqRecord<true> fastqIstring = fastqString;

    REQUIRE(fastqIstring.name == "TestSequence");
    REQUIRE(Codec::to_string(fastqIstring.seq) == "ACGT");
    REQUIRE(fastqIstring.qual == "!@?#");
  }

  SECTION("Parsing from Stream") {
    std::istringstream stream("@TestSequence\nACGT\n+\n!@?#");
    FastqRecord<> fastqString;
    stream >> fastqString;

    REQUIRE(fastqString.name == "TestSequence");
    REQUIRE(fastqString.seq == "ACGT");
    REQUIRE(fastqString.qual == "!@?#");
  }

  SECTION("Writing to Stream") {
    FastqRecord<> fastqString;
    fastqString.name = "TestSequence";
    fastqString.seq = "ACGT";
    fastqString.qual = "!@?#";

    std::ostringstream stream;
    stream << fastqString;

    std::string expectedOutput = "@TestSequence\nACGT\n+\n!@?#";
    REQUIRE(stream.str() == expectedOutput);
  }

  SECTION("Implicit Conversion from istring to string") {
    FastqRecord<true> fastqIstring;
    fastqIstring.name = "TestSequence";
    fastqIstring.seq = "0123"_s;
    fastqIstring.qual = "!@?#";

    // Test implicit conversion from istring to string
    FastqRecord<> fastqString = fastqIstring;

    REQUIRE(fastqString.name == "TestSequence");
    REQUIRE(fastqString.seq == "ACGT");
    REQUIRE(fastqString.qual == "!@?#");
  }

  SECTION("Parsing from Stream with Different Encodings") {
    std::istringstream stringStream("@TestSequence\nACGT\n+\n!@?#");
    FastqRecord<> fastqString;
    stringStream >> fastqString;

    std::istringstream istringStream("@TestSequence\nACGT\n+\n!@?#");
    FastqRecord<true> fastqIstring;
    istringStream >> fastqIstring;

    REQUIRE(fastqString.name == fastqIstring.name);
    REQUIRE(fastqString.seq == Codec::to_string(fastqIstring.seq));
    REQUIRE(fastqString.qual == fastqIstring.qual);
  }

  SECTION("Writing to Stream with Different Encodings") {
    FastqRecord<> fastqString;
    fastqString.name = "TestSequence";
    fastqString.seq = "ACGT";
    fastqString.qual = "!@?#";

    FastqRecord<true> fastqIstring;
    fastqIstring.name = "TestSequence";
    fastqIstring.seq = "0123"_s;
    fastqIstring.qual = "!@?#";

    std::ostringstream stringStream;
    stringStream << fastqString;

    std::ostringstream istringStream;
    istringStream << fastqIstring;

    std::string expectedStringOutput = "@TestSequence\nACGT\n+\n!@?#";
    std::string expectedIstringOutput = "@TestSequence\nACGT\n+\n!@?#";

    REQUIRE(stringStream.str() == expectedStringOutput);
    REQUIRE(istringStream.str() == expectedIstringOutput);
  }
}

TEST_CASE("Fastq basic I/O") {
  SECTION("Single record") {
    FastqRecord record;

    std::istringstream iss{R"(@SRR001666.1
GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9IC
)"};

    std::ostringstream oss;

    while (iss >> record) oss << record << '\n';

    REQUIRE(record.name == "SRR001666.1");
    REQUIRE(record.seq == "GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACC");
    REQUIRE(record.qual == "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9IC");
  }

  SECTION("Multiple records") {
    FastqRecord record1, record2, record3, record4;

    std::istringstream iss{R"(@alpha1 comm
AAAAAAAAA
+
IIIIIIIII
@alpha2 comm
TTTTTTTTT
+
IIIIIIIII
@alpha3 comm
CCCCCCCCC
+
IIIIIIIII
@alpha4 comm
GGGGGGGGG
+
IIIIIIIII


)"};

    iss >> record1 >> record2 >> record3 >> record4;

    REQUIRE(record1.name == "alpha1");
    REQUIRE(record1.seq == "AAAAAAAAA");
    REQUIRE(record1.qual == "IIIIIIIII");

    REQUIRE(record2.name == "alpha2");
    REQUIRE(record2.seq == "TTTTTTTTT");
    REQUIRE(record2.qual == "IIIIIIIII");

    REQUIRE(record3.name == "alpha3");
    REQUIRE(record3.seq == "CCCCCCCCC");
    REQUIRE(record3.qual == "IIIIIIIII");

    REQUIRE(record4.name == "alpha4");
    REQUIRE(record4.seq == "GGGGGGGGG");
    REQUIRE(record4.qual == "IIIIIIIII");
  }

  SECTION("SHORTQUALITY.FASTQ") {
    std::istringstream iss{R"(@SEQ1
ACfTACGTACGTAGCTGATCGATCGTACGTAGCTGACA
+
SHORTQUALITY:)
@SEQ2
NNNNNCGTACGTAGCTGATCGATCGTACGTAGCTGACA
+
!!!!!AIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII)"};

    FastqRecord record1, record2;
    while (iss >> record1 >> record2)
      ;
    REQUIRE(record1.name == "SEQ1");
    REQUIRE(record1.seq == "ACfTACGTACGTAGCTGATCGATCGTACGTAGCTGACA");
    REQUIRE(record1.qual == "SHORTQUALITY:)");

    REQUIRE(record2.name == "SEQ2");
    REQUIRE(record2.seq == "NNNNNCGTACGTAGCTGATCGATCGTACGTAGCTGACA");
    REQUIRE(record2.qual == "!!!!!AIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII");
  }

  SECTION("OLD_SANGER.FASTQ") {
    std::istringstream iss{R"(@SANGER_FASTQ
ACGTGCTAGCTAGCTGATCGTACGTAGCTGACT
ACGTGCTAGCTAGCTGATCGTACGTAGCTGACT
ACGTGCTAGCTAGCTGATCGTACGTAGCTGACT
ACGTGCTAGCTAGCTGATCGTACGTAGCTGACT
+
999999999999999897989999999989889
999664999999999897989999999989889
999999199999999897989999999989889
999999911999999897989999999989889)"};

    FastqRecord record;

    while (iss >> record)
      ;

    REQUIRE(record.name == "SANGER_FASTQ");
    REQUIRE(
        record.seq ==
        "ACGTGCTAGCTAGCTGATCGTACGTAGCTGACTACGTGCTAGCTAGCTGATCGTACGTAGCTGACTACGT"
        "GCTAGCTAGCTGATCGTACGTAGCTGACTACGTGCTAGCTAGCTGATCGTACGTAGCTGACT");
    REQUIRE(
        record.qual ==
        "9999999999999998979899999999898899996649999999998979899999999898899999"
        "99199999999897989999999989889999999911999999897989999999989889");
  }
}


TEST_CASE("FASTQ FILE I/O") {
  SECTION("READFILE1.FASTQ") {
    std::string_view expected_result{R"(@A00709:43:HYG25DSXX:1:1101:3640:1000
GCATTCACCCTGGTCGGGTCGGCGTTGTAATCTGCCTGGACCAGACTACGCACTGTCGGTGGGGTGGCGGCGCGGGAAACGTCATGTCGC
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
@A00709:43:HYG25DSXX:1:1101:6189:1000
TCGGAACCTCGTCCACGATTTGCGGAGCCGCGTTCGCGACCAGGCGGTCCTTGCCCACCAACTGCAGGGTCATCAAGTAGCCCCCGGGGC
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFF:FFFFFFFFFFFFFFFFFFFFF,FFFFFFF
@A00709:43:HYG25DSXX:1:1101:10818:1000
CCAACAATGCTTACGTTCACACCCAACGCCCGAACCCTATGACGGTAGGCAAGTTAAGGCGGGCTTTTTTGCGGATTTACGTAAAGCGGC
+
FFFFFFFFF:FFFFFFFFFFFFFFFFFF:FFFF:F:FF,FFFFFFFFFFFFFF:FFFFFFFFFFFFF:FFF:FF::F:FFF,F:FFFFFF
@A00709:43:HYG25DSXX:1:1101:14705:1000
GTTGTGCCGGTAATAACATTTGTTATTGAGAGGGCCCTCCTGCGATTGGCTTGATGGTTCCGTAAGTGTGAAATGTCACTCCGTTATCGA
+
FFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFF
@A00709:43:HYG25DSXX:1:1101:14814:1000
ATCGTAGGCAAAAGTTCCAACAAGATCTGTTGTCCCGCCGTCTGTCCCAAGGGGATTAACTTCACCCCCGCCCCAATTAAATTGCTCGCC
+
FFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
@A00709:43:HYG25DSXX:1:1101:15935:1000
CTCTTCCTTGTCGGCGGTGGAGAAGCAGAGGCAGAAAAGGTCGCGCTCGTAGGGGATGCCGTCCTCGATCGCCAGACGCTCGCTGGCGCG
+
FFFFF,FFFFFFFFFFFFFFFF:FFFFF:FFFFFFFFF,FFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFF
@A00709:43:HYG25DSXX:1:1101:17960:1000
CCCGCACGCGGATGTTCGGATGGTTCCGCACGGCCTCGCTGAGGGCGGTCTGGATCACCAATCCGGTGGCGTCCTTGGCGTGGTAGATGC
+
FF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:
)"};
    CHECK_FASTQ_IDENTITY(data_path / "test1.fastq", expected_result);
  }

  SECTION("READFILE2.FASTQ") {
    std::string_view expected_result{R"(@SEQ1
ACGTACGTACGTAGCTGATCGATCGTACGTAGCTGACA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@SEQ2
NNNNNCGTACGTAGCTGATCGATCGTACGTAGCTGACA
+
!!!!!AIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@SEQ3
ACGTACGTACGTAGCTGATCGATCGTACGTAGCTGACN
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
)"};
    CHECK_FASTQ_IDENTITY(data_path / "test2.fastq", expected_result);
  }

  SECTION("READFILE3.FASTQ") {
    std::string_view expected_result{R"(@A00709:43:HYG25DSXX:1:1101:3640:1000
GCATTCACCCTGGTCGGGTCGGCGTTGTAATCTGCCTGGACCAGACTACGCACTGTCGGTGGGGTGGCGGCGCGGGAAACGTCATGTCGC
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
@A00709:43:HYG25DSXX:1:1101:3640:1000
ACCACGACTTACGTGGATGGCAATGTGACGGTCGGAACCGAATACGAATATCGCGTGGAGCGCACGGGGTCGTCCTTCGACGGAAATGCC
+
FFF:FFFFFFFF:FFFFF:FFFFFFFFFFFF:FFFFFFFFFFFFFF:FFF,FFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFF:FFFFF
@A00709:43:HYG25DSXX:1:1101:6189:1000
TCGGAACCTCGTCCACGATTTGCGGAGCCGCGTTCGCGACCAGGCGGTCCTTGCCCACCAACTGCAGGGTCATCAAGTAGCCCCCGGGGC
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFF:FFFFFFFFFFFFFFFFFFFFF,FFFFFFF
@A00709:43:HYG25DSXX:1:1101:6189:1000
TCGACGCGGTCGCTAAGTTCACAGATGCCGTTCAGATGGACATCGCCGTGCGGTTGTCCCTTACGCCCGACGACCCAGGTGCGGTCTCGC
+
FFFFFFFF:FFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFF:FFFFFFFFFFFFFFF:F:FFFFFFFFFFFFFFFFFFF:FFFFFFF
@A00709:43:HYG25DSXX:1:1101:10818:1000
CCAACAATGCTTACGTTCACACCCAACGCCCGAACCCTATGACGGTAGGCAAGTTAAGGCGGGCTTTTTTGCGGATTTACGTAAAGCGGC
+
FFFFFFFFF:FFFFFFFFFFFFFFFFFF:FFFF:F:FF,FFFFFFFFFFFFFF:FFFFFFFFFFFFF:FFF:FF::F:FFF,F:FFFFFF
@A00709:43:HYG25DSXX:1:1101:10818:1000
CGCCAAAGGCGCGCCGCAGCGACCGATGCACAGCGGCCCGCGCTCCGACGCGCCGCCCAGAGGCCCACGGGGACCGCGCGGCCCGGGCAA
+
FFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFF
@A00709:43:HYG25DSXX:1:1101:14705:1000
GTTGTGCCGGTAATAACATTTGTTATTGAGAGGGCCCTCCTGCGATTGGCTTGATGGTTCCGTAAGTGTGAAATGTCACTCCGTTATCGA
+
FFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFF
@A00709:43:HYG25DSXX:1:1101:14705:1000
ATGCCATAGATCATTGAGATTTCAAGGTTGGAAGGAGAGAAGTATATATGTTAATACCACGAAGAAATCTGGTGAAATTTGGTTGGGTTA
+
,FF,FFFF:F:FFFFFFFFFFFFF:FFFF:F:FF,FFFFFFFFFFF:FFFF,FFFFFFFFF,F::FFFFFFFFF,,F,F::,FFFFFFFF
@A00709:43:HYG25DSXX:1:1101:14814:1000
ATCGTAGGCAAAAGTTCCAACAAGATCTGTTGTCCCGCCGTCTGTCCCAAGGGGATTAACTTCACCCCCGCCCCAATTAAATTGCTCGCC
+
FFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
@A00709:43:HYG25DSXX:1:1101:14814:1000
CCTGTGCACAACTCCCTTTAACTAACCCCAGAATCATTATTAAGACTTCAACAACTAGCAAGTCCTATCTTGCCTGTCGGGAGTTACAGC
+
F:FFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFF:F:FFFFFFFFFFFFFFFFFFFFFFF:FFFFFF
@A00709:43:HYG25DSXX:1:1101:15935:1000
CTCTTCCTTGTCGGCGGTGGAGAAGCAGAGGCAGAAAAGGTCGCGCTCGTAGGGGATGCCGTCCTCGATCGCCAGACGCTCGCTGGCGCG
+
FFFFF,FFFFFFFFFFFFFFFF:FFFFF:FFFFFFFFF,FFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFF
@A00709:43:HYG25DSXX:1:1101:15935:1000
GACGCCGTTCGACCAGCGCCACGCCATGCGCAGCCCGCGGATCTTCGACGTCGTGGCGTCGTTCCCCAAGCCGGTGATCGCCATGATCAA
+
:FFFFFFFFFFFFFFFFFFFFFFFF:FFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF,FFFFFFF:FFFFFFFFFFFF
@A00709:43:HYG25DSXX:1:1101:17960:1000
CCCGCACGCGGATGTTCGGATGGTTCCGCACGGCCTCGCTGAGGGCGGTCTGGATCACCAATCCGGTGGCGTCCTTGGCGTGGTAGATGC
+
FF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:
@A00709:43:HYG25DSXX:1:1101:17960:1000
ATCATCGGCCTGGCGCCGCCGGAGGAAGGGGATTCGCCCGAGCTGCTGGCCGCGGATATCGAAGCCGCGGGCGCGGGCCTGTGCCGGACC
+
FFFF:FFFFFFFFFFFFFFFFF,FF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFF:FFFFFFF
)"};
    CHECK_FASTQ_IDENTITY(data_path / "test3.fastq", expected_result);
  }
  SECTION("READFILE4.FASTQ") {
    std::string_view expected_result{R"()"};
    CHECK_FASTQ_IDENTITY(data_path / "test4.fastq", expected_result);
  }
  SECTION("READFILE5.FASTQ") {
    std::string_view expected_result{
      R"(@M03029:193:000000000-D2B6C:1:1101:15779:1330
CTTAGAAG
+
>AABBBDF
@M03029:193:000000000-D2B6C:1:1101:15821:1338
CTTAGAAG
+
ABCCCFFF
@M03029:193:000000000-D2B6C:1:1101:15389:1341
CTTAGAAG
+
>AABBFFF
@M03029:193:000000000-D2B6C:1:1101:15215:1344
CTTAGAAG
+
>AABBFBB
@M03029:193:000000000-D2B6C:1:1101:15519:1345
CTTAGAAG
+
>ABCCFFF
@M03029:193:000000000-D2B6C:1:1101:16240:1357
CTTAGAAG
+
3AAAADFF
@M03029:193:000000000-D2B6C:1:1101:15715:1360
CTTAGAAG
+
3AAAAF5B
@M03029:193:000000000-D2B6C:1:1101:15337:1360
CTTAGAAG
+
>1>A1131
@M03029:193:000000000-D2B6C:1:1101:15736:1362
CTTAGAAG
+
BBBCCFFF
@M03029:193:000000000-D2B6C:1:1101:15887:1362
CTTAGAAG
+
ABCBBFFF
@M03029:193:000000000-D2B6C:1:1101:16608:1365
CTTAGAAG
+
AABCBFFF
@M03029:193:000000000-D2B6C:1:1101:15106:1373
CTTAGAAG
+
>ABBBFFF
@M03029:193:000000000-D2B6C:1:1101:15404:1374
CTTAGAAG
+
3>AAABFF
@M03029:193:000000000-D2B6C:1:1101:15382:1375
CTTAGAAG
+
AABCCFFF
@M03029:193:000000000-D2B6C:1:1101:16235:1375
CTTAGAAG
+
>ABBAFFF
@M03029:193:000000000-D2B6C:1:1101:14878:1378
CTTAGAAG
+
1>>A11B1
@M03029:193:000000000-D2B6C:1:1101:14777:1378
CTTAGAAG
+
>ABCCFFF
@M03029:193:000000000-D2B6C:1:1101:16596:1379
CTTAGAAG
+
ABBBBFFF
@M03029:193:000000000-D2B6C:1:1101:16832:1379
CTTAAAAG
+
AABBAFFF
@M03029:193:000000000-D2B6C:1:1101:16907:1383
CTTAGAAG
+
AAABCFFF
)"};
    CHECK_FASTQ_IDENTITY(data_path / "test5.fastq", expected_result);
  }

  SECTION("READFILE6.FASTQ") {
    std::string_view expected_result{R"(@SEQ1
A
+
I
@SEQ2
N
+
!
@SEQ3
A
+
I
)"};
    CHECK_FASTQ_IDENTITY(data_path / "test6.fastq", expected_result);
  }

  SECTION("READFILE7.FASTQ") {
    std::string_view expected_result{R"(@SRR001666.1
GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCGGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCGGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCGGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCGGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9IC
)"};
    CHECK_FASTQ_IDENTITY(data_path / "test7.fastq", expected_result);
  }

  SECTION("READFILE8.FASTQ") {
    std::string_view expected_result{R"(@SEQ_ID
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
)"};
    CHECK_FASTQ_IDENTITY(data_path / "test8.fastq", expected_result);
  }

  SECTION("READFILE9.FASTQ") {
    std::string_view expected_result{R"(@MN00537:51:000H2K25G:1:11101:2213:1092
CTCCAGTCCTTACTCCCATATCTAACCTCTTACCCCTACNTCATAGGTANACATTTTAATGAAT
+
FFFFFFFFFFFFAFFFFFFFF=FFFFAFFFFFFF/AFFF#FFFFFFFFF#FFFFFFFF
@MN00537:51:000H2K25G:1:11101:2213:1092
CTCCAGTCCTTACTCCCATATCTAACCTCTTACCCCTACNTCATAGGTANACATTTTAATGAAT
+
FFFFFFFFFFFFAFFFFFFFF=FFFFAFFFFFFF/AFFF#FFFFFFFFF#FFFFFFFF
@MN00537:51:000H2K25G:1:11101:2213:1092
CTCCAGTCCTTACTCCCATATCTAACCTCTTACCCCTACNTCATAGGTANACATTTTAATGAAT
+
FFFFFFFFFFFFAFFFFFFFF=FFFFAFFFFFFF/AFFF#FFFFFFFFF#FFFFFFFF
@SRR001666.1
GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCGGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCGGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCGGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCGGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9IC
)"};
    CHECK_FASTQ_IDENTITY(data_path / "test9.fastq", expected_result);
  }

  SECTION("READFILE10.FASTQ") {
    std::string_view expected_result{R"(@empSEQ

+

)"};
    CHECK_FASTQ_IDENTITY(data_path / "test10.fastq", expected_result);
  }
}
