#pragma once
#include <biovoltron/file_io/core/record.hpp>
#include <biovoltron/utility/interval.hpp>
#include <boost/algorithm/string/iter_find.hpp>
#include <boost/algorithm/string/finder.hpp>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <filesystem>
#include <fstream>


namespace biovoltron {


/**
 * @ingroup file_io
 * @brief A format represent six column bed file, where the last two columns are gene_type and gene_name for annotations
 */
struct BedSixRecord : Record {
  std::string seqid;
  std::uint32_t start{};
  std::uint32_t end{};
  char strand{};
  /**
   * @brief gene_type, such as protein_coding, tRNA, or miRNA.
   */
  std::string gene_type;
  /**
   * @brief gene_name, such as MNX1 (protein_coding), miR101-a-1-3p (miRNA) 
   */
  std::string gene_name;

  auto
  operator<=>(const BedSixRecord& other) const noexcept {
    return std::tie(seqid, start, end)
           <=> std::tie(other.seqid, other.start, other.end);
  }

  operator auto() const { return Interval{seqid, start, end, strand}; }
};

namespace bedsixreader {
  using bedsix = BedSixRecord;
  using bedsix_v = std::vector<bedsix>;

  template<class Str>
  std::vector<std::string> split(Str&& str, const std::string& delimiter )
  {
    std::vector<std::string> results;
    boost::iter_split(results, str, boost::algorithm::first_finder(delimiter));
    return results;
  }


  template<class Path>
  void read_mirtrondb_gff(Path&& input_file_path, bedsix_v& container) {
    /// Please download mirtronDB GFF file for your species from http://mirtrondb.cp.utfpr.edu.br/download.php
    /**
     * 0         1      2       3       4   5       6       7       8
     * seq_id	source	type	start	end	score	strand	phase	atributtes
     * 7 	 mirtronDB 	 SO:0001034 	 66381669 	 66381690 	 . 	 - 	 . 	 mirtronDB_id:mmu-mirtron-1855-3p;organism:mmusculus;prec_mature:m;arm_steem:3p;miRBase_id:mmu-mir-7057;
     * 7 	 mirtronDB 	 SO:0001034 	 66381702 	 66381719 	 . 	 - 	 . 	 mirtronDB_id:mmu-mirtron-1855-5p;organism:mmusculus;prec_mature:m;arm_steem:5p;miRBase_id:mmu-mir-7057;
     */

    if (!std::filesystem::exists(input_file_path)) {
      std::string file_name = input_file_path;
      throw std::runtime_error("read_mirtrondb_gff file " + file_name + " does not exist!!");
    }
    std::ifstream file(input_file_path);
    std::string line;
    std::string tips_msg = R"(
Start reading the mirtrondb gff annotation file.
If any `ERROR` or `Segmentation fault` occurred, 
the reason might be the breaking change of the file format of the annotation.
Please check the following URL for mirtrondb gff file format:
http://mirtrondb.cp.utfpr.edu.br/
    )";
    std::cout << tips_msg << std::endl;

    while(std::getline(file, line)) {
      // std::cout << line << std::endl;
      if (line.at(0) == '#') continue;
      std::vector<std::string> content = split(line, "\t");
      if (content[0] == "seq_id") continue;
      std::for_each(content.begin(), content.end(), [](auto&& str){
        boost::trim_right(str);
        boost::trim_left(str);
      });
      std::string arm = split(split(content[8], "arm_steem:")[1], ";")[0];
      std::string name = split(split(content[8], "miRBase_id:")[1], ";")[0];
      if (name == "-") continue;
      name = name + "-" + arm;
      // std::vector<std::string> names = split(attributes[1], ";");
      std::string seqid = std::string{"chr"} + content[0];
      // std::string name = names[0];
      container.emplace_back( 
        bedsix {
          {}, // for Record
          seqid,
          static_cast<uint32_t>(std::stoul(content[3])),
          static_cast<uint32_t>(std::stoul(content[4])),
          *content[6].begin(),
          "mirtron",
          name
        }
      );
    }
    file.close();
  }

  template<class Path>
  bedsix_v read_mirtrondb_gff_into(Path&& input_file_path) {
    bedsix_v results;
    read_mirtrondb_gff(input_file_path, results);
    return results;
  }

  template<class Path>
  void read_mirbase_gff(Path&& input_file_path, bedsix_v& container) {
    /// Please download mirBase gff file for your species from https://www.mirbase.org/ftp.shtml
    /// e.g. for mm10, wget https://www.mirbase.org/ftp/CURRENT/genomes/mmu.gff3
    /**
     * chr1	.	miRNA	12426016	12426038	.	+	.	ID=MIMAT0025084;Alias=MIMAT0025084;Name=mmu-miR-6341;Derives_from=MI0021869
     */
    if (!std::filesystem::exists(input_file_path)) {
      std::string file_name = input_file_path;
      throw std::runtime_error("read_mirbase_gff file " + file_name + " does not exist!!");
    }
    std::ifstream file(input_file_path);
    std::string line;
    std::string tips_msg = R"(
Start reading the mirbase gff annotation file.
If any `ERROR` or `Segmentation fault` occurred, 
the reason might be the breaking change of the file format of the annotation.
Please check the following URL for mirbase gff file format:
https://www.mirbase.org/
    )";
    std::cout << tips_msg << std::endl;

    while(std::getline(file, line)) {
      if (line.at(0) == '#') continue;
      std::vector<std::string> content = split(line, "\t");
      if (content[2] != "miRNA") continue;
      std::vector<std::string> attribute = split(content[8], "Name=");
      std::string name = split(attribute[1], ";")[0];
      
      container.emplace_back( 
        bedsix {
          {}, // for Record
          content[0],
          static_cast<uint32_t>(std::stoul(content[3])),
          static_cast<uint32_t>(std::stoul(content[4])),
          *content[6].begin(),
          content[2],
          name
        }
      );
    }
    file.close();
  }

  template<class Path>
  bedsix_v read_mirbase_gff_into(Path&& input_file_path) {
    bedsix_v results;
    read_mirbase_gff(input_file_path, results);
    return results;
  }
  
  template<class Path>
  void read_gencode_gtf(Path&& input_file_path, bedsix_v& container) {
    /// gtf example from gencode
    /**
     * chr1	HAVANA	exon	3073253	3074322	.	+	.	gene_id "ENSMUSG00000102693.1"; transcript_id "ENSMUST00000193812.1"; gene_type "TEC"; gene_name "4933401J01Rik"; transcript_type "TEC"; transcript_name "4933401J01Rik-201"; exon_number 1; exon_id "ENSMUSE00001343744.1"; level 2; transcript_support_level "NA"; tag "basic"; havana_gene "OTTMUSG00000049935.1"; havana_transcript "OTTMUST00000127109.1";
     * chr1	ENSEMBL	gene	3102016	3102125	.	+	.	gene_id "ENSMUSG00000064842.1"; gene_type "snRNA"; gene_name "Gm26206"; level 3;
     * chr1	ENSEMBL	transcript	3102016	3102125	.	+	.	gene_id "ENSMUSG00000064842.1"; transcript_id "ENSMUST00000082908.1"; gene_type "snRNA"; gene_name "Gm26206"; transcript_type "snRNA"; transcript_name "Gm26206-201"; level 3; transcript_support_level "NA"; tag "basic";
     */
    if (!std::filesystem::exists(input_file_path)) {
      std::string file_name = input_file_path;
      throw std::runtime_error("read_gencode_gtf file " + file_name + " does not exist!!");
    }
    std::ifstream file(input_file_path);
    std::string line;
    std::string tips_msg = R"(
Start reading the gencode gtf annotation file.
If any `ERROR` or `Segmentation fault` occurred, 
the reason might be the breaking change of the file format of the annotation.
Please check the following URL for gencode gtf file format:
https://www.gencodegenes.org/pages/data_format.html
    )";

    std::cout << tips_msg << std::endl;
    while(std::getline(file, line)) {
      if (line.at(0) == '#') continue;
      std::vector<std::string> content = split(line, "\t");
      std::vector<std::string> split_gene_type = split(content[8], "gene_type \"");
      std::vector<std::string> split_gene_name = split(content[8], "gene_name \"");
      std::string gene_type = split(split_gene_type[1], "\";")[0];
      std::string gene_name = split(split_gene_name[1], "\";")[0];
      container.emplace_back( 
        bedsix {
          {}, // for Record
          content[0],
          static_cast<uint32_t>(std::stoul(content[3])),
          static_cast<uint32_t>(std::stoul(content[4])),
          *content[6].begin(),
          gene_type,
          gene_name
        }
      );
    }
    file.close();
  }

  template<class Path>
  bedsix_v read_gencode_gtf_into(Path&& input_file_path) {
    bedsix_v results;
    read_gencode_gtf(input_file_path, results);
    return results;
  }
} // namespace bedsixreader

}  // namespace biovoltron