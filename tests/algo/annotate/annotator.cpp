#include <biovoltron/algo/annotate/annotator.hpp>
#include <biovoltron/file_io/sam.hpp>
#include <biovoltron/file_io/gff.hpp>
#include <biovoltron/file_io/bed.hpp>
#include <biovoltron/file_io/vcf.hpp>
#include <fstream>
#include <filesystem>
#include <catch.hpp>

using namespace biovoltron;

const auto data_path = std::filesystem::path{DATA_PATH};

struct Feature : Interval {
  std::string name;
  std::string type;
  std::uint32_t id;
};

bool overlap(Interval a, Interval b) {
  return a.chrom == b.chrom && a.begin < b.end && a.end > b.begin;
}

TEST_CASE("Annotator") {
  SECTION("Insert an object at a location") {
    auto genes = Annotator<std::string>{};
    genes.insert_at("gene1", Interval{"chr1:5-15"});
    genes.insert_at("gene2", Interval{"chr1:2-10"});
    genes.insert_at("gene3", Interval{"chr1:20-30"});
    genes.insert_at("gene4", Interval{"-chr2:2-10"});
    genes.index();

    // result returns a vector of object, sorted by begin
    auto results = genes.find(Interval{"chr1:6-12"});
    REQUIRE(results.size() == 2);
    REQUIRE(results[0] == "gene2");
    REQUIRE(results[1] == "gene1");

    results = genes.find(Interval{"-chr2:6-9"});
    REQUIRE(results.size() == 1);
    REQUIRE(results[0] == "gene4");

    results = genes.find(Interval{"chr2:6-9"});
    REQUIRE(results.empty());

    results = genes.find(Interval{"chr2:0-2"});
    REQUIRE(results.empty());

    results = genes.find(Interval{"chr999:5-6"});
    REQUIRE(results.empty());
  }

  SECTION("Insert an object that has location info") {
    auto genes = Annotator<Feature>{};
    auto gene1 = Feature{"chr1:5-15", "gene1", "gene", 5566};
    auto gene2 = Feature{"chr1:2-10", "gene2", "gene", 5567};
    auto gene3 = Feature{"chr1:20-30", "gene3", "gene", 5568};
    auto gene4 = Feature{"-chr2:2-10", "gene4", "gene", 5569};
    genes.insert(gene1);
    genes.insert(gene2);
    genes.insert(gene3);
    genes.insert(gene4);
    genes.index();

    auto results = genes.find(Interval{"chr1:6-12"});
    REQUIRE(results.size() == 2);
    REQUIRE(std::ranges::is_sorted(results));
    REQUIRE(results[0] == gene2);
    REQUIRE(results[1] == gene1);

    results = genes.find(Interval{"-chr2:6-9"});
    REQUIRE(results.size() == 1);
    REQUIRE(results[0] == gene4);

    results = genes.find(Interval{"+chr2:6-9"});
    REQUIRE(results.empty());

    results = genes.find(Interval{"chr999:5-6"});
    REQUIRE(results.empty());
  }

  SECTION("Throw when find() before index()") {
    auto genes = Annotator<std::string>{};
    genes.insert_at("gene1", Interval{"chr1:2-10"});
    REQUIRE_THROWS(genes.find(Interval{"chr1:4-6"}));
  }

  SECTION("WARNING") {
    auto genes = Annotator<int>{};
    genes.insert_at(1, "gene"); // XXX: this is equivalent to genes.insert_at(1, "gene:0-4294967295")
  }

  SECTION("GFF and SAM") {
    std::fstream fin(data_path/"gene.gff");
    GffRecord gff_record;
    Annotator<GffRecord> features;
    while (fin >> gff_record)
      features.insert(gff_record);
    features.index();

    fin.close();
    fin.clear();
    fin.open(data_path/"test2.sam");
    SamRecord sam_record;
    while (fin >> sam_record) {
      auto results = features.find(sam_record);
      REQUIRE_FALSE(results.empty());
      for (const auto& result : results) 
        REQUIRE(Interval{sam_record}.overlaps(result));
    }
  }

  SECTION("BED and SAM") {
    std::fstream fin(data_path/"gene.bed");
    BedRecord bed_record;
    Annotator<BedRecord> features;
    while (fin >> bed_record)
      features.insert(bed_record);
    features.index();
    fin.close();
    fin.clear();
    fin.open(data_path/"test2.sam");
    SamRecord sam_record;
    while (fin >> sam_record) {
      auto results = features.find(sam_record);
      REQUIRE_FALSE(results.empty());
      for (const auto& result : results)
        REQUIRE(Interval{sam_record}.overlaps(result));
    }
  }

  SECTION("VCF and SAM") {
    std::fstream fin(data_path/"test.vcf");
    VcfRecord vcf_record;
    VcfHeader vcf_header;
    Annotator<VcfRecord> features;
    fin >> vcf_header;
    while (fin >> vcf_record) 
      features.insert(vcf_record);
    features.index();

    fin.close();
    fin.clear();
    fin.open(data_path/"test2.sam");
    SamRecord sam_record;
    while (fin >> sam_record) {
      auto results = features.find(sam_record);
      for (const auto& result : results) {
        REQUIRE(Interval{sam_record}.overlaps(result));
      }
    }
  }
}
