#include <biovoltron/file_io/gff.hpp>
#include <catch.hpp>
#include <sstream>
#include <iostream>

using namespace biovoltron;

TEST_CASE("gff") {
  GffRecord r;
  std::istringstream iss{
      "ctg123\t.\tmRNA\t10000\t15000\t0\t+\t0\tID=mrna0002;Parent=operon001;Name=subsonicsquirrel"};
  iss >> r;
  REQUIRE(r.seqid == "ctg123");
  REQUIRE(r.source == ".");
  REQUIRE(r.type == "mRNA");
  REQUIRE(r.start == 10000);
  REQUIRE(r.end == 15000);
  REQUIRE(r.score == 0);
  REQUIRE(r.strand == '+');
  REQUIRE(r.phase == 0);
  REQUIRE(r.attrs == "ID=mrna0002;Parent=operon001;Name=subsonicsquirrel");
  REQUIRE(Interval{r} == Interval{"ctg123", 10000 - 1, 15000, '+'});
}

TEST_CASE("gff missing fields") {
  GffRecord r;
  std::istringstream iss{
      ".\t.\t.\t.\t.\t.\t.\t.\t."};
  iss >> r;
  REQUIRE(r.seqid == ".");
  REQUIRE(r.source == ".");
  REQUIRE(r.type == ".");
  REQUIRE(r.start == 0);
  REQUIRE(r.end == 0);
  REQUIRE(r.score == 0);
  REQUIRE(r.strand == '.');
  REQUIRE(r.phase == 0);
  REQUIRE(r.attrs == ".");
}

TEST_CASE("gff Equal to Comparison") {
  GffRecord rec1,rec2;

  std::istringstream iss1{"ctg123\t.\tmRNA\t10000\t15000\t.\t+\t.\tID=mrna0002;Parent=operon001;Name=subsonicsquirrel"};
  std::istringstream iss2{"ctg123\t.\tmRNA\t10000\t15000\t.\t+\t.\tID=mrna0002;Parent=operon001;Name=subsonicsquirrel"};
  iss1 >> rec1;
  iss2 >> rec2;
  REQUIRE(false==(rec1<rec2));
  REQUIRE(true==(rec1==rec2));
  REQUIRE(false==(rec1>rec2));



}

TEST_CASE("gff Less than Comparison") {
  GffRecord rec1,rec2;

  /* Seqid */
  std::istringstream iss1{"btg123\t.\tmRNA\t10000\t15000\t.\t+\t.\tID=mrna0002;Parent=operon001;Name=subsonicsquirrel"};
  std::istringstream iss2{"ctg123\t.\tmRNA\t10000\t15000\t.\t+\t.\tID=mrna0002;Parent=operon001;Name=subsonicsquirrel"};
  iss1 >> rec1;
  iss2 >> rec2;
  REQUIRE(true==(rec1<rec2));
  REQUIRE(false==(rec1==rec2));
  REQUIRE(false==(rec1>rec2));
  iss1.clear();
  iss2.clear();

  /* Start */
  iss1.str("ctg123\t.\tmRNA\t9999\t15000\t.\t+\t.\tID=mrna0002;Parent=operon001;Name=subsonicsquirrel");
  iss2.str("ctg123\t.\tmRNA\t10000\t15000\t.\t+\t.\tID=mrna0002;Parent=operon001;Name=subsonicsquirrel");
  iss1 >> rec1;
  iss2 >> rec2;
  REQUIRE(true==(rec1<rec2));
  REQUIRE(false==(rec1==rec2));
  REQUIRE(false==(rec1>rec2));
  iss1.clear();
  iss2.clear();

  /* End */
  iss1.str("ctg123\t.\tmRNA\t10000\t14999\t.\t+\t.\tID=mrna0002;Parent=operon001;Name=subsonicsquirrel");
  iss2.str("ctg123\t.\tmRNA\t10000\t15000\t.\t+\t.\tID=mrna0002;Parent=operon001;Name=subsonicsquirrel");
  iss1 >> rec1;
  iss2 >> rec2;
  REQUIRE(true==(rec1<rec2));
  REQUIRE(false==(rec1==rec2));
  REQUIRE(false==(rec1>rec2));
}

TEST_CASE("gff Greater than Comparison") {
  GffRecord rec1,rec2;

  /* Seqid */
  std::istringstream iss1{"dtg123\t.\tmRNA\t10000\t15000\t.\t+\t.\tID=mrna0002;Parent=operon001;Name=subsonicsquirrel"};
  std::istringstream iss2{"ctg123\t.\tmRNA\t10000\t15000\t.\t+\t.\tID=mrna0002;Parent=operon001;Name=subsonicsquirrel"};
  iss1 >> rec1;
  iss2 >> rec2;
  REQUIRE(false==(rec1<rec2));
  REQUIRE(false==(rec1==rec2));
  REQUIRE(true==(rec1>rec2));
  iss1.clear();
  iss2.clear();

  /* Start */
  iss1.str("ctg123\t.\tmRNA\t10001\t15000\t.\t+\t.\tID=mrna0002;Parent=operon001;Name=subsonicsquirrel");
  iss2.str("ctg123\t.\tmRNA\t10000\t15000\t.\t+\t.\tID=mrna0002;Parent=operon001;Name=subsonicsquirrel");
  iss1 >> rec1;
  iss2 >> rec2;
  REQUIRE(false==(rec1<rec2));
  REQUIRE(false==(rec1==rec2));
  REQUIRE(true==(rec1>rec2));
  iss1.clear();
  iss2.clear();

  /* End */
  iss1.str("ctg123\t.\tmRNA\t10000\t15001\t.\t+\t.\tID=mrna0002;Parent=operon001;Name=subsonicsquirrel");
  iss2.str("ctg123\t.\tmRNA\t10000\t15000\t.\t+\t.\tID=mrna0002;Parent=operon001;Name=subsonicsquirrel");
  iss1 >> rec1;
  iss2 >> rec2;
  REQUIRE(false==(rec1<rec2));
  REQUIRE(false==(rec1==rec2));
  REQUIRE(true==(rec1>rec2));
  iss1.clear();
  iss2.clear();

}
