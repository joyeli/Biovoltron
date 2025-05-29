#include <biovoltron/file_io/bam.hpp>
#include <vector>
#include <sstream>
#include <spdlog/spdlog.h>
#include <catch.hpp>

using namespace biovoltron;

std::filesystem::path data_path = DATA_PATH;
const auto in = data_path / "test.bam";
auto bam_index = std::filesystem::path(in.string() + ".bai");
auto out = data_path / "a.bam";

template<bool Encoded = false>
auto same(SamRecord<Encoded>& r, std::string_view expected) {
  std::stringstream ss;
  REQUIRE(true);
  ss << r;
  return ss.str() == expected;
}

template<bool Encoded = false>
auto test_encoded() {
  IBamStream fin(in);
  SamHeader h1, h2;
  SamRecord<Encoded> s1, s2, s3;
  SECTION("Basic input") {
    fin >> h1;
    fin >> s1;
    REQUIRE(h1.lines.size() == 9);
    REQUIRE(same(s1, "HWI-ST486:305:C0RH5ACXX:1:2104:8917:83075	99	1	1150	40	43S13M6872N45M	=	8048	170	CAGACAGGAACTAGCAATGCTTGAAATCAAGAACTTGAATTGAAATAGTTTTTTACTGGATCAGAGACTACTCAATATCCCCAAACTTGGAAATTAGTTTG	CCCFFFFFHHHHHJIJJJJIJJJJJJJJJJJIJJJJJJIJJJJJJJJJBGGJJJJJJJJIIJJJJJJJIHHHHHFFFFFFFDEDEDDDDDCDCCDDDCCED	MD:Z:2G55	NH:i:1	HI:i:1	NM:i:1	SM:i:40	XQ:i:40	X2:i:0	XS:A:-	"));
    while (fin >> s3);
    REQUIRE(fin.eof());
    SECTION("to_begin()") {
      REQUIRE(same(s1, "HWI-ST486:305:C0RH5ACXX:1:2104:8917:83075	99	1	1150	40	43S13M6872N45M	=	8048	170	CAGACAGGAACTAGCAATGCTTGAAATCAAGAACTTGAATTGAAATAGTTTTTTACTGGATCAGAGACTACTCAATATCCCCAAACTTGGAAATTAGTTTG	CCCFFFFFHHHHHJIJJJJIJJJJJJJJJJJIJJJJJJIJJJJJJJJJBGGJJJJJJJJIIJJJJJJJIHHHHHFFFFFFFDEDEDDDDDCDCCDDDCCED	MD:Z:2G55	NH:i:1	HI:i:1	NM:i:1	SM:i:40	XQ:i:40	X2:i:0	XS:A:-	"));
      fin.to_begin();
      fin >> s2;
      REQUIRE(s1 == s2);
    }
    SECTION("Set region") {
      REQUIRE(fin.set_region("1", 0, 20000));
      fin >> s2;
      REQUIRE(s1 == s2);
    }
  }
  
  SECTION("Load index file") {
    REQUIRE(fin.load_index(bam_index));
    REQUIRE(fin.is_indexed());
    REQUIRE(fin.on_sequential());
  }
  SECTION("Set unmapped") {
    REQUIRE(fin.set_unmapped());
    REQUIRE(fin.on_unmapped());
    // TODO: need a bam that contain unmapped read to test read property
  }
  
  SECTION("close") {
    fin.close();
    // TODO: if there is any input after close the stream, the failbit
    //       should be set
    // fin >> s4;
    // assert(!fin.good());
  }
}

TEST_CASE("IBamStream") {
  SECTION("Test encoded SamRecord") {
    test_encoded<true>();
  }
  SECTION("Test no encoded SamRecord") {
    test_encoded<false>();
  }
}

TEST_CASE("OBamStream") {
  SECTION("conditional") {
    OBamStream fout(out);
    REQUIRE(fout.is_open());
  }
  SECTION("Basic output") {
    SamHeader h1, h2;
    SamRecord<false> s1, s2, s3;
    auto foo = [&](bool gen_idx = false) {  
      {
        IBamStream fin(in);
        fin >> h1;
        fin >> s1;
      }
      {
        OBamStream fout(out, gen_idx);
        fout << h1;
        fout << s1;
      }
      {
        IBamStream fin(out);
        fin >> h2;
        fin >> s2;
        REQUIRE(h1 == h2);
        REQUIRE(s1 == s2);
      }
    };
    foo(false);
    foo(true);
  }
}