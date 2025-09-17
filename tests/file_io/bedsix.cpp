#include <biovoltron/file_io/bedsix.hpp>
#include <catch.hpp>
#include <filesystem>
#include <sstream>
#include <iostream>

using namespace biovoltron;

const auto root_path = std::filesystem::path{DATA_PATH};
const auto data_path = root_path / "tailor";
const auto mirbase_path = data_path / "mmu_mirBase_v22.gff3";
const auto mirtron_path = data_path / "mmu_mirtron.gff";
const auto main_path = data_path / "mmu_gencode_main.gtf";
const auto pseudo_path = data_path / "mmu_gencode_pseudos.gtf";
const auto polyA_path = data_path / "mmu_gencode_polyAs.gtf";
const auto tRNAs_path = data_path / "mmu_gencode_tRNAs.gtf";

TEST_CASE("BedSixRecord::Parsing - Parses BED6 records", "[BedSixRecord]") {
  BedSixRecord r;
  std::istringstream iss{
      "chr7\t12\t127472363\t+\tmiRNA\tmiR92a-1-3p"};
  iss >> r;
  REQUIRE(r.seqid == "chr7");
  REQUIRE(r.start == 12);
  REQUIRE(r.end == 127472363);
  REQUIRE(r.strand == '+');
  REQUIRE(r.gene_type == "miRNA");
  REQUIRE(r.gene_name == "miR92a-1-3p");
  REQUIRE(Interval{r} == Interval{"chr7", 12, 127472363, '+'});
}

TEST_CASE("BedSixReader::Reading - Reads BED6 files", "[BedSixRecord]") {
  using namespace bedsixreader;
  SECTION("read from mirbase gff") {
    std::vector<BedSixRecord> records;
    INFO("Reading mirbase annotation file failed.");
    INFO("The reason might be the breaking change of the file format of the annotation.");
    INFO("For more information, please checkout: ");
    INFO("https://www.mirbase.org");
    read_mirbase_gff(mirbase_path, records);

    auto first = records.front();
    REQUIRE(first.seqid == "chr1");
    REQUIRE(first.start == 12426015);
    REQUIRE(first.end == 12426038);
    REQUIRE(first.strand == '+');
    REQUIRE(first.gene_type == "miRNA");
    REQUIRE(first.gene_name == "mmu-miR-6341");
    REQUIRE(Interval{first} == Interval{"chr1", 12426015, 12426038, '+'});

    auto last = records.back();
    REQUIRE(last.seqid == "chr1");
    REQUIRE(last.start == 133726618);
    REQUIRE(last.end == 133726640);
    REQUIRE(last.strand == '-');
    REQUIRE(last.gene_type == "miRNA");
    REQUIRE(last.gene_name == "mmu-miR-6903-5p");
    REQUIRE(Interval{last} == Interval{"chr1", 133726618, 133726640, '-'});
  }
  SECTION("read from mirbase gff into") {
    std::vector<BedSixRecord> records = read_mirbase_gff_into(mirbase_path);

    auto first = records.front();
    REQUIRE(first.seqid == "chr1");
    REQUIRE(first.start == 12426015);
    REQUIRE(first.end == 12426038);
    REQUIRE(first.strand == '+');
    REQUIRE(first.gene_type == "miRNA");
    REQUIRE(first.gene_name == "mmu-miR-6341");
    REQUIRE(Interval{first} == Interval{"chr1", 12426015, 12426038, '+'});
  }
  SECTION("read from mirtronDB gff") {
    std::vector<BedSixRecord> records;
    INFO("Reading mirtronDB annotation file failed.");
    INFO("The reason might be the breaking change of the file format of the annotation.");
    INFO("For more information, please checkout: ");
    INFO("http://mirtrondb.cp.utfpr.edu.br/");
    read_mirtrondb_gff(mirtron_path, records);

    auto first = records.front();
    REQUIRE(first.seqid == "chr7");
    REQUIRE(first.start == 66381668);
    REQUIRE(first.end == 66381690);
    REQUIRE(first.strand == '-');
    REQUIRE(first.gene_type == "mirtron");
    REQUIRE(first.gene_name == "mmu-mir-7057-3p");
    REQUIRE(Interval{first} == Interval{"chr7", 66381668, 66381690, '-'});

    auto last = records.back();
    REQUIRE(last.seqid == "chr8");
    REQUIRE(last.start == 71631046);
    REQUIRE(last.end == 71631067);
    REQUIRE(last.strand == '-');
    REQUIRE(last.gene_type == "mirtron");
    REQUIRE(last.gene_name == "mmu-mir-6769b-3p");
    REQUIRE(Interval{last} == Interval{"chr8", 71631046, 71631067, '-'});

    
  }
  SECTION("read from mirtronDB gff into") {
    std::vector<BedSixRecord> records = read_mirtrondb_gff_into(mirtron_path);

    auto first = records.front();
    REQUIRE(first.seqid == "chr7");
    REQUIRE(first.start == 66381668);
    REQUIRE(first.end == 66381690);
    REQUIRE(first.strand == '-');
    REQUIRE(first.gene_type == "mirtron");
    REQUIRE(first.gene_name == "mmu-mir-7057-3p");
    REQUIRE(Interval{first} == Interval{"chr7", 66381668, 66381690, '-'});
  }
  SECTION("read from gencode gtf - main") {
    std::vector<BedSixRecord> records;
    INFO("Reading gencode gtf annotation file failed.");
    INFO("The reason might be the breaking change of the file format of the annotation.");
    INFO("For more information, please checkout: ");
    INFO("https://www.gencodegenes.org/pages/data_format.html");

    read_gencode_gtf(main_path, records);
    auto first = records.front();
    REQUIRE(first.seqid == "chr1");
    REQUIRE(first.start == 3073252);
    REQUIRE(first.end == 3074322);
    REQUIRE(first.strand == '+');
    REQUIRE(first.gene_type == "TEC");
    REQUIRE(first.gene_name == "4933401J01Rik");
    REQUIRE(Interval{first} == Interval{"chr1", 3073252, 3074322, '+'});

    auto last = records.back();
    REQUIRE(last.seqid == "chr2");
    REQUIRE(last.start == 58567297);
    REQUIRE(last.end == 58792971);
    REQUIRE(last.strand == '+');
    REQUIRE(last.gene_type == "protein_coding");
    REQUIRE(last.gene_name == "Upp2");
    REQUIRE(Interval{last} == Interval{"chr2", 58567297, 58792971, '+'});
  }
  SECTION("read from gencode gtf into") {
    std::vector<BedSixRecord> records = read_gencode_gtf_into(main_path);

    auto last = records.back();
    REQUIRE(last.seqid == "chr2");
    REQUIRE(last.start == 58567297);
    REQUIRE(last.end == 58792971);
    REQUIRE(last.strand == '+');
    REQUIRE(last.gene_type == "protein_coding");
    REQUIRE(last.gene_name == "Upp2");
    REQUIRE(Interval{last} == Interval{"chr2", 58567297, 58792971, '+'});
  }
  SECTION("read from gencode gtf - polyAs") {
    std::vector<BedSixRecord> records;
    INFO("Reading gencode gtf annotation file failed.");
    INFO("The reason might be the breaking change of the file format of the annotation.");
    INFO("For more information, please checkout: ");
    INFO("https://www.gencodegenes.org/pages/data_format.html");

    read_gencode_gtf(polyA_path, records, "polyA_site");
    auto first = records.front();
    REQUIRE(first.seqid == "chr1");
    REQUIRE(first.start == 3214480);
    REQUIRE(first.end == 3214482);
    REQUIRE(first.strand == '-');
    REQUIRE(first.gene_type == "polyA_site");
    REQUIRE(first.gene_name == "744347");
    REQUIRE(Interval{first} == Interval{"chr1", 3214480, 3214482, '-'});

    auto last = records.back();
    REQUIRE(last.seqid == "chr5");
    REQUIRE(last.start == 30659729);
    REQUIRE(last.end == 30659731);
    REQUIRE(last.strand == '+');
    REQUIRE(last.gene_type == "polyA_site");
    REQUIRE(last.gene_name == "703508");
    REQUIRE(Interval{last} == Interval{"chr5", 30659729, 30659731, '+'});
  }
  SECTION("read from gencode gtf - pseudos") {
    std::vector<BedSixRecord> records;
    INFO("Reading gencode gtf annotation file failed.");
    INFO("The reason might be the breaking change of the file format of the annotation.");
    INFO("For more information, please checkout: ");
    INFO("https://www.gencodegenes.org/pages/data_format.html");

    read_gencode_gtf(pseudo_path, records, "transcript");
    auto first = records.front();
    REQUIRE(first.seqid == "chr1");
    REQUIRE(first.start == 3252733);
    REQUIRE(first.end == 3253236);
    REQUIRE(first.strand == '+');
    REQUIRE(first.gene_type == "pseudogene");
    REQUIRE(first.gene_name == "PGOMOU00000274813");
    REQUIRE(Interval{first} == Interval{"chr1", 3252733, 3253236, '+'});

    auto last = records.back();
    REQUIRE(last.seqid == "chr3");
    REQUIRE(last.start == 94604410);
    REQUIRE(last.end == 94604784);
    REQUIRE(last.strand == '-');
    REQUIRE(last.gene_type == "pseudogene");
    REQUIRE(last.gene_name == "PGOMOU00000276566");
    REQUIRE(Interval{last} == Interval{"chr3", 94604410, 94604784, '-'});
  }
  SECTION("read from gencode gtf - tRNAs") {
    std::vector<BedSixRecord> records;
    INFO("Reading gencode gtf annotation file failed.");
    INFO("The reason might be the breaking change of the file format of the annotation.");
    INFO("For more information, please checkout: ");
    INFO("https://www.gencodegenes.org/pages/data_format.html");

    read_gencode_gtf(tRNAs_path, records, "tRNA");
    auto first = records.front();
    REQUIRE(first.seqid == "chr1");
    REQUIRE(first.start == 112349388);
    REQUIRE(first.end == 112349461);
    REQUIRE(first.strand == '+');
    REQUIRE(first.gene_type == "Pseudo_tRNA");
    REQUIRE(first.gene_name == "NULL");
    REQUIRE(Interval{first} == Interval{"chr1", 112349388, 112349461, '+'});

    auto last = records.back();
    REQUIRE(last.seqid == "chr11");
    REQUIRE(last.start == 113168793);
    REQUIRE(last.end == 113168861);
    REQUIRE(last.strand == '-');
    REQUIRE(last.gene_type == "Ala_tRNA");
    REQUIRE(last.gene_name == "NULL");
    REQUIRE(Interval{last} == Interval{"chr11", 113168793, 113168861, '-'});
  }
}

TEST_CASE("BedSixReader::ErrorHandling - Handles file reading errors", "[BedSixRecord]") {
  using namespace bedsixreader;
  SECTION("file does not exist") {
    std::vector<BedSixRecord> records;
    auto exception_msg_gencode = [](auto&& str){
      return "read_gencode_gtf file " + std::string{str} + " does not exist!!";
    };
    auto exception_msg_mirtrondb = [](auto&& str){
      return "read_mirtrondb_gff file " + std::string{str} + " does not exist!!";
    };
    auto exception_msg_mirbase = [](auto&& str){
      return "read_mirbase_gff file " + std::string{str} + " does not exist!!";
    };
    auto nonexist_path = data_path / "nonexist_file";
    REQUIRE_THROWS_WITH(read_gencode_gtf(nonexist_path, records), exception_msg_gencode(nonexist_path));
    REQUIRE_THROWS_WITH(read_mirtrondb_gff(nonexist_path, records), exception_msg_mirtrondb(nonexist_path));
    REQUIRE_THROWS_WITH(read_mirbase_gff(nonexist_path, records), exception_msg_mirbase(nonexist_path));
  }
}