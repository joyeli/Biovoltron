#include <biovoltron/file_io/bed.hpp>
#include <catch.hpp>
#include <sstream>
#include <iostream>

using namespace biovoltron;

TEST_CASE("bed") {
  BedRecord r;
  std::istringstream iss{
      "chr7\t12\t127472363\tPos1\t0\t+\t127471196\t127472363\t255,0,0\t3\t354,109,1189\t0,739,1347,"};
  iss >> r;
  REQUIRE(r.chrom == "chr7");
  REQUIRE(r.start == 12);
  REQUIRE(r.end == 127472363);
  REQUIRE(r.name == "Pos1");
  REQUIRE(r.score == 0);
  REQUIRE(r.strand == '+');
  REQUIRE(r.thick_start == 127471196);
  REQUIRE(r.thick_end == 127472363);
  REQUIRE(r.item_rgb == "255,0,0");
  REQUIRE(r.block_count == 3);
  REQUIRE(r.block_sizes == "354,109,1189");
  REQUIRE(r.block_starts == "0,739,1347,");
  REQUIRE(Interval{r} == Interval{"chr7", 12, 127472363, '+'});
}

TEST_CASE("BED missing fields"){
  BedRecord r;
  std::istringstream iss{
    "chr7\t127471196\t127472363\tPos1"};
  iss >> r;
  REQUIRE(r.chrom == "chr7");
  REQUIRE(r.start == 127471196);
  REQUIRE(r.end == 127472363);
  REQUIRE(r.name == "Pos1");
  REQUIRE(r.score == 0);
  REQUIRE(r.strand == 0);
  REQUIRE(r.thick_start == 0);
  REQUIRE(r.thick_end == 0);
  REQUIRE(r.item_rgb == "0,0,0");
  REQUIRE(r.block_count == 0);
  REQUIRE(r.block_sizes == "0");
  REQUIRE(r.block_starts == "0");
}

TEST_CASE("BED Equal to Comparison"){
  BedRecord rec1,rec2;
  
  /* Equal */
  std::istringstream iss1{
    "chr7\t12\t127472363\tPos1\t0\t+\t127471196\t127472363\t255,0,0\t3\t354,109,1189\t0,739,1347,"
    };
  std::istringstream iss2{
    "chr7\t12\t127472363\tPos1\t0\t+\t127471196\t127472363\t255,0,0\t3\t354,109,1189\t0,739,1347,"
    };
  iss1>>rec1;
  iss2>>rec2;
  
  REQUIRE(false==(rec1<rec2));
  REQUIRE(false==(rec1>rec2));
  REQUIRE(true==(rec1==rec2));
}

TEST_CASE("BED Less than Comparison"){
  BedRecord rec1,rec2;

  /* Chrom */
  std::istringstream iss1{
    "chr1\t12\t127472362\tPos1\t0\t+\t127471196\t127472363\t255,0,0\t3\t354,109,1189\t0,739,1347,"
    };
  std::istringstream iss2{
    "chr2\t12\t127472362\tPos1\t0\t+\t127471196\t127472363\t255,0,0\t3\t354,109,1189\t0,739,1347,"
    };
  iss1>>rec1;
  iss2>>rec2;
  
  REQUIRE(true==(rec1<rec2));
  REQUIRE(false==(rec1>rec2));
  REQUIRE(false==(rec1==rec2));

  iss1.clear();
  iss2.clear();

  iss1.str("chr1\t12\t127472362\tPos1\t0\t+\t127471196\t127472363\t255,0,0\t3\t354,109,1189\t0,739,1347,");
  iss2.str("chr11\t12\t127472362\tPos1\t0\t+\t127471196\t127472363\t255,0,0\t3\t354,109,1189\t0,739,1347,");

  iss1>>rec1;
  iss2>>rec2;
  
  REQUIRE(true==(rec1<rec2));
  REQUIRE(false==(rec1>rec2));
  REQUIRE(false==(rec1==rec2));

  iss1.clear();
  iss2.clear();

  iss1.str("chr1\t12\t127472362\tPos1\t0\t+\t127471196\t127472363\t255,0,0\t3\t354,109,1189\t0,739,1347,");
  iss2.str("chrX\t12\t127472362\tPos1\t0\t+\t127471196\t127472363\t255,0,0\t3\t354,109,1189\t0,739,1347,");

  iss1>>rec1;
  iss2>>rec2;
  
  REQUIRE(true==(rec1<rec2));
  REQUIRE(false==(rec1>rec2));
  REQUIRE(false==(rec1==rec2));

  iss1.clear();
  iss2.clear();

  iss1.str("chrX\t12\t127472362\tPos1\t0\t+\t127471196\t127472363\t255,0,0\t3\t354,109,1189\t0,739,1347,");
  iss2.str("chrY\t12\t127472362\tPos1\t0\t+\t127471196\t127472363\t255,0,0\t3\t354,109,1189\t0,739,1347,");

  iss1>>rec1;
  iss2>>rec2;
  
  REQUIRE(true==(rec1<rec2));
  REQUIRE(false==(rec1>rec2));
  REQUIRE(false==(rec1==rec2));

  iss1.clear();
  iss2.clear();

  /* Start */
  iss1.str("chr1\t11\t127472362\tPos1\t0\t+\t127471196\t127472363\t255,0,0\t3\t354,109,1189\t0,739,1347,");
  iss2.str("chr1\t12\t127472362\tPos1\t0\t+\t127471196\t127472363\t255,0,0\t3\t354,109,1189\t0,739,1347,");
  iss1>>rec1;
  iss2>>rec2;
  
  REQUIRE(true==(rec1<rec2));
  REQUIRE(false==(rec1>rec2));
  REQUIRE(false==(rec1==rec2));

  iss1.clear();
  iss2.clear();

  /* End */

  iss1.str("chr1\t12\t127472362\tPos1\t0\t+\t127471196\t127472363\t255,0,0\t3\t354,109,1189\t0,739,1347,");
  iss2.str("chr1\t12\t127472363\tPos1\t0\t+\t127471196\t127472363\t255,0,0\t3\t354,109,1189\t0,739,1347,");
  iss1>>rec1;
  iss2>>rec2;
  
  REQUIRE(true==(rec1<rec2));
  REQUIRE(false==(rec1>rec2));
  REQUIRE(false==(rec1==rec2));

  iss1.clear();
  iss2.clear();

}

TEST_CASE("BED Greater than Comparison"){
  BedRecord rec1,rec2;

  /* Chrom */
  std::istringstream iss1{
    "chr2\t12\t127472362\tPos1\t0\t+\t127471196\t127472363\t255,0,0\t3\t354,109,1189\t0,739,1347,"
    };
  std::istringstream iss2{
    "chr1\t12\t127472362\tPos1\t0\t+\t127471196\t127472363\t255,0,0\t3\t354,109,1189\t0,739,1347,"
    };
  iss1>>rec1;
  iss2>>rec2;
  
  REQUIRE(false==(rec1<rec2));
  REQUIRE(true==(rec1>rec2));
  REQUIRE(false==(rec1==rec2));

  iss1.clear();
  iss2.clear();

  iss1.str("chr11\t12\t127472362\tPos1\t0\t+\t127471196\t127472363\t255,0,0\t3\t354,109,1189\t0,739,1347,");
  iss2.str("chr1\t12\t127472362\tPos1\t0\t+\t127471196\t127472363\t255,0,0\t3\t354,109,1189\t0,739,1347,");

  iss1>>rec1;
  iss2>>rec2;
  
  REQUIRE(false==(rec1<rec2));
  REQUIRE(true==(rec1>rec2));
  REQUIRE(false==(rec1==rec2));

  iss1.clear();
  iss2.clear();

  iss1.str("chrX\t12\t127472362\tPos1\t0\t+\t127471196\t127472363\t255,0,0\t3\t354,109,1189\t0,739,1347,");
  iss2.str("chr1\t12\t127472362\tPos1\t0\t+\t127471196\t127472363\t255,0,0\t3\t354,109,1189\t0,739,1347,");

  iss1>>rec1;
  iss2>>rec2;
  
  REQUIRE(false==(rec1<rec2));
  REQUIRE(true==(rec1>rec2));
  REQUIRE(false==(rec1==rec2));

  iss1.clear();
  iss2.clear();

  iss1.str("chrY\t12\t127472362\tPos1\t0\t+\t127471196\t127472363\t255,0,0\t3\t354,109,1189\t0,739,1347,");
  iss2.str("chrX\t12\t127472362\tPos1\t0\t+\t127471196\t127472363\t255,0,0\t3\t354,109,1189\t0,739,1347,");

  iss1>>rec1;
  iss2>>rec2;
  
  REQUIRE(false==(rec1<rec2));
  REQUIRE(true==(rec1>rec2));
  REQUIRE(false==(rec1==rec2));

  iss1.clear();
  iss2.clear();

  /* Start */
  iss1.str("chr1\t12\t127472362\tPos1\t0\t+\t127471196\t127472363\t255,0,0\t3\t354,109,1189\t0,739,1347,");
  iss2.str("chr1\t11\t127472362\tPos1\t0\t+\t127471196\t127472363\t255,0,0\t3\t354,109,1189\t0,739,1347,");
  iss1>>rec1;
  iss2>>rec2;
  
  REQUIRE(false==(rec1<rec2));
  REQUIRE(true==(rec1>rec2));
  REQUIRE(false==(rec1==rec2));

  iss1.clear();
  iss2.clear();

  /* End */

  iss1.str("chr1\t12\t127472363\tPos1\t0\t+\t127471196\t127472363\t255,0,0\t3\t354,109,1189\t0,739,1347,");
  iss2.str("chr1\t12\t127472362\tPos1\t0\t+\t127471196\t127472363\t255,0,0\t3\t354,109,1189\t0,739,1347,");
  iss1>>rec1;
  iss2>>rec2;

  REQUIRE(false==(rec1<rec2));
  REQUIRE(true==(rec1>rec2));
  REQUIRE(false==(rec1==rec2));

  iss1.clear();
  iss2.clear();

}

TEST_CASE("BED Header")
{
  BedHeader Bh;
  std::istringstream iss{
      "browser position chr7:127471196-127495720\n"
      "browser hide all\n"
      "track name=HbVar type=bedDetail description=\"HbVar custom track\" db=hg19 visibility=3 url=\"http://globin.bx.psu.edu/cgi-bin/hbvar/query_vars3?display_format=page&mode=output&id=$$\"\n"};  
      
  iss >> Bh;

  std::ostringstream oss;
  auto p_cout_streambuf = std::cout.rdbuf();
  std::cout.rdbuf(oss.rdbuf());

  std::cout<<Bh<<std::endl;
  std::cout.rdbuf(p_cout_streambuf); // restore
  REQUIRE(oss.str() == iss.str());

}

TEST_CASE("bed_graph") {
  BedGraphRecord r;
  std::istringstream iss{"chr19\t49302000\t49302300\t-1.0"};
  iss >> r;
  REQUIRE(r.chrom == "chr19");
  REQUIRE(r.start == 49302000);
  REQUIRE(r.end == 49302300);
  REQUIRE(r.score == -1);
}

TEST_CASE("BedGraph missing fields") {
  BedGraphRecord r;
  std::istringstream iss{".\t.\t.\t."};
  iss >> r;
  REQUIRE(r.chrom == ".");
  REQUIRE(r.start == 0);
  REQUIRE(r.end == 0);
  REQUIRE(r.score == 0);
}

TEST_CASE("BedGraph Header")
{
  BedHeader Bh;
  std::istringstream iss{
  "browser position chr19:49302001-49304701\n"
  "browser hide all\n"
  "browser pack refGene encodeRegions\n"
  "browser full altGraph\n"
  "#	300 base wide bar graph, autoScale is on by default == graphing\n"
  "#	limits will dynamically change to always show full range of data\n"
  "#	in viewing window, priority = 20 positions this as the second graph\n"
  "#	Note, zero-relative, half-open coordinate system in use for bedGraph format\n"
  "track type=bedGraph name=\"BedGraph Format\" description=\"BedGraph format\" visibility=full color=200,100,0 altColor=0,100,200 priority=20\n"};  
      
  iss >> Bh;

  std::ostringstream oss;
  auto p_cout_streambuf = std::cout.rdbuf();
  std::cout.rdbuf(oss.rdbuf());

  std::cout<<Bh<<std::endl;

  std::cout.rdbuf(p_cout_streambuf); // restore
  REQUIRE(oss.str() == iss.str());

}
