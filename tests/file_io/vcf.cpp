#include <biovoltron/file_io/vcf.hpp>
#include <catch.hpp>
#include <sstream>
#include <filesystem>

using namespace biovoltron;
using namespace std::string_literals;

template<typename T>
using Iter = std::istream_iterator<T>;

const auto data_path = std::filesystem::path{DATA_PATH};

TEST_CASE("vcf") {
  VcfRecord r;
  std::istringstream iss{
      "20\t1110696\trs6040355\tA\tG,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4"};
  iss >> r;
  REQUIRE(r.chrom == "20");
  REQUIRE(r.pos == 1110696);
  REQUIRE(r.id == "rs6040355");
  REQUIRE(r.ref == "A");
  REQUIRE(r.alt == "G,T");
  REQUIRE(r.qual == 67);
  REQUIRE(r.filter == "PASS");
  REQUIRE(r.info == "NS=2;DP=10;AF=0.333,0.667;AA=T;DB");
  REQUIRE(r.format == "GT:GQ:DP:HQ");
  REQUIRE(r.samples.size() == 3);
  REQUIRE(r.samples[0] == "1|2:21:6:23,27");
  REQUIRE(r.samples[1] == "2|1:2:0:18,2");
  REQUIRE(r.samples[2] == "2/2:35:4");  
  REQUIRE(Interval{r} == Interval{"20", 1110696 - 1, 1110696, '+'});
}
