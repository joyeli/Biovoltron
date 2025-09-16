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
    // std::cout << tips_msg << std::endl;

    while(std::getline(file, line)) {
      // std::cout << line << std::endl;
      if (line.empty()) continue;
      if (line.at(0) == '#') continue;
      std::vector<std::string> content = split(line, "\t");
      if (content[0] == "seq_id" ||
          content[3].empty() || content[4].empty() ||
          content[3] == "  " || content[4] == "  ") continue;
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
          static_cast<uint32_t>(std::stoul(content[3])) - 1, // 1-based (gff) to 0-based (bed)
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
    // std::cout << tips_msg << std::endl;

    while(std::getline(file, line)) {
      if(line.empty()) continue;
      if (line.at(0) == '#') continue;
      std::vector<std::string> content = split(line, "\t");
      if (content[2] != "miRNA") continue;
      std::vector<std::string> attribute = split(content[8], "Name=");
      std::string name = split(attribute[1], ";")[0];

      // FIXME: This will break the naming convention of miRNA, don't do this
      // // append ID tag (like '3' of 'MIMAT0041637_3') to name
      // attribute = split(content[8], "ID=");
      // std::string id = split(attribute[1], ";")[0];
      // if(auto underscore_pos = id.find('_'); underscore_pos != std::string::npos){
      //   name += id.substr(underscore_pos);
      // }
      
      container.emplace_back( 
        bedsix {
          {}, // for Record
          content[0],
          static_cast<uint32_t>(std::stoul(content[3])) - 1, // 1-based (gff) to 0-based (bed)
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
  void read_gencode_gtf(Path&& input_file_path, bedsix_v& container, const std::string& feature_type = "gene") {
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

    // std::cout << tips_msg << std::endl;
    while(std::getline(file, line)) {
      if(line.empty()) continue;
      if (line.at(0) == '#') continue;
      std::vector<std::string> content = split(line, "\t");
      if (content[2] != feature_type) continue;
      std::vector<std::string> split_gene_type = split(content[8], "gene_type \"");
      std::vector<std::string> split_gene_name = split(content[8], "gene_name \"");
      std::string gene_type = split(split_gene_type[1], "\";")[0];
      std::string gene_name = split(split_gene_name[1], "\";")[0];
      container.emplace_back( 
        bedsix {
          {}, // for Record
          content[0],
          static_cast<uint32_t>(std::stoul(content[3])) - 1, // 1-based (gff) to 0-based (bed)
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

  template<class Path>
  void read_gencode_gff(Path&& input_file_path, bedsix_v& container) {
    /// gff example from gencode
    /**
     * chr1	HAVANA	gene	11869	14409	.	+	.	ID=ENSG00000290825.1;gene_id=ENSG00000290825.1;gene_type=lncRNA;gene_name=DDX11L2;level=2;tag=overlaps_pseudogene
     * chr1	HAVANA	transcript	11869	14409	.	+	.	ID=ENST00000456328.2;Parent=ENSG00000290825.1;gene_id=ENSG00000290825.1;transcript_id=ENST00000456328.2;gene_type=lncRNA;gene_name=DDX11L2;transcript_type=lncRNA;transcript_name=DDX11L2-202;level=2;transcript_support_level=1;tag=basic,Ensembl_canonical;havana_transcript=OTTHUMT00000362751.1
     * chr1	HAVANA	exon	11869	12227	.	+	.	ID=exon:ENST00000456328.2:1;Parent=ENST00000456328.2;gene_id=ENSG00000290825.1;transcript_id=ENST00000456328.2;gene_type=lncRNA;gene_name=DDX11L2;transcript_type=lncRNA;transcript_name=DDX11L2-202;exon_number=1;exon_id=ENSE00002234944.1;level=2;transcript_support_level=1;tag=basic,Ensembl_canonical;havana_transcript=OTTHUMT00000362751.1
     */
    if (!std::filesystem::exists(input_file_path)) {
      std::string file_name = input_file_path;
      throw std::runtime_error("read_gencode_gff file " + file_name + " does not exist!!");
    }
    std::ifstream file(input_file_path);
    std::string line;
    std::string tips_msg = R"(
Start reading the gencode gff annotation file.
If any `ERROR` or `Segmentation fault` occurred, 
the reason might be the breaking change of the file format of the annotation.
Please check the following URL for gencode gff file format:
https://www.gencodegenes.org/pages/data_format.html
)";

    // std::cout << tips_msg << std::endl;
    while(std::getline(file, line)) {
      if(line.empty()) continue;
      if (line.at(0) == '#') continue;
      std::vector<std::string> content = split(line, "\t");
      if (content[2] != "gene") continue;
      std::vector<std::string> split_gene_type = split(content[8], "gene_type=");
      std::vector<std::string> split_gene_name = split(content[8], "gene_name=");
      std::string gene_type = split(split_gene_type[1], ";")[0];
      std::string gene_name = split(split_gene_name[1], ";")[0];
      container.emplace_back( 
        bedsix {
          {}, // for Record
          content[0],
          static_cast<uint32_t>(std::stoul(content[3])) - 1, // 1-based (gff) to 0-based (bed)
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
  bedsix_v read_gencode_gff_into(Path&& input_file_path) {
    bedsix_v results;
    read_gencode_gff(input_file_path, results);
    return results;
  }
  
  template<class Path>
  void read_gencode_trna_gtf(Path&& input_file_path, bedsix_v& container) {
    /// gtf example from gencode
    /**
     * chr1	ENSEMBL	tRNA	7930279	7930348	.	-	.	gene_id "33323"; transcript_id "33323"; gene_type "Pseudo_tRNA"; gene_name "33323"; transcript_type "Pseudo_tRNA"; transcript_name "33323"; level 3;
     * chr1	ENSEMBL	tRNA	16520585	16520658	.	-	.	gene_id "20658"; transcript_id "20658"; gene_type "Asn_tRNA"; gene_name "20658"; transcript_type "Asn_tRNA"; transcript_name "20658"; level 3;
     * chr1	ENSEMBL	tRNA	16532398	16532471	.	-	.	gene_id "20657"; transcript_id "20657"; gene_type "Asn_tRNA"; gene_name "20657"; transcript_type "Asn_tRNA"; transcript_name "20657"; level 3;
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

    // std::cout << tips_msg << std::endl;
    while(std::getline(file, line)) {
      if(line.empty()) continue;
      if (line.at(0) == '#') continue;
      std::vector<std::string> content = split(line, "\t");
      std::vector<std::string> split_gene_type = split(content[8], "gene_type \"");
      std::vector<std::string> split_gene_name = split(content[8], "transcript_id \"");
      std::string gene_type = split(split_gene_type[1], "_")[0];
      std::string gene_name = split(split_gene_name[1], "\";")[0];
      std::string strand = content[6];
      std::size_t start = std::stoul(content[3]) - 1; // 1-based (gff) to 0-based (bed)
      std::size_t end = std::stoul(content[4]);
      std::size_t mid = (end - start) / 2 + start;
      gene_name = gene_type + "-" + gene_name; // Asn_20657
      container.emplace_back( 
        bedsix {
          {}, // for Record
          content[0],
          static_cast<uint32_t>(start),
          static_cast<uint32_t>(mid),
          *content[6].begin(),
          "tRF",
          gene_name + (strand == "+" ? "-5p" : "-3p")
        }
      );
      container.emplace_back( 
        bedsix {
          {}, // for Record
          content[0],
          static_cast<uint32_t>(mid),
          static_cast<uint32_t>(end),
          *content[6].begin(),
          "tRF",
          gene_name + (strand == "+" ? "-3p" : "-5p")
        }
      );
    }
    file.close();
  }

  template<class Path>
  bedsix_v read_gencode_trna_gtf_into(Path&& input_file_path) {
    bedsix_v results;
    read_gencode_trna_gtf(input_file_path, results);
    return results;
  }

  template<class Path>
  void read_gencode_trna_gff(Path&& input_file_path, bedsix_v& container) {
    /// gff example from gencode
    /**
     * chr1	ENSEMBL	tRNA	7930279	7930348	.	-	.	ID=33323;gene_id=33323;transcript_id=33323;gene_type=Pseudo_tRNA;gene_name=33323;transcript_type=Pseudo_tRNA;transcript_name=33323;level=3
     * chr1	ENSEMBL	tRNA	16520585	16520658	.	-	.	ID=20658;gene_id=20658;transcript_id=20658;gene_type=Asn_tRNA;gene_name=20658;transcript_type=Asn_tRNA;transcript_name=20658;level=3
     * chr1	ENSEMBL	tRNA	16532398	16532471	.	-	.	ID=20657;gene_id=20657;transcript_id=20657;gene_type=Asn_tRNA;gene_name=20657;transcript_type=Asn_tRNA;transcript_name=20657;level=3
     */
    if (!std::filesystem::exists(input_file_path)) {
      std::string file_name = input_file_path;
      throw std::runtime_error("read_gencode_gff file " + file_name + " does not exist!!");
    }
    std::ifstream file(input_file_path);
    std::string line;
    std::string tips_msg = R"(
Start reading the gencode gff annotation file.
If any `ERROR` or `Segmentation fault` occurred, 
the reason might be the breaking change of the file format of the annotation.
Please check the following URL for gencode gff file format:
https://www.gencodegenes.org/pages/data_format.html
)";

    // std::cout << tips_msg << std::endl;
    while(std::getline(file, line)) {
      if(line.empty()) continue;
      if (line.at(0) == '#') continue;
      std::vector<std::string> content = split(line, "\t");
      std::vector<std::string> split_gene_type = split(content[8], "gene_type=");
      std::vector<std::string> split_gene_name = split(content[8], "transcript_id=");
      std::string gene_type = split(split_gene_type[1], "_")[0];
      std::string gene_name = split(split_gene_name[1], ";")[0];
      std::string strand = content[6];
      std::size_t start = std::stoul(content[3]) - 1; // 1-based (gff) to 0-based (bed)
      std::size_t end = std::stoul(content[4]);
      std::size_t mid = (end - start) / 2 + start;
      gene_name = gene_type + "-" + gene_name; // Asn_20657
      container.emplace_back( 
        bedsix {
          {}, // for Record
          content[0],
          static_cast<uint32_t>(start),
          static_cast<uint32_t>(mid),
          *content[6].begin(),
          "tRF",
          gene_name + (strand == "+" ? "-5p" : "-3p")
        }
      );
      container.emplace_back( 
        bedsix {
          {}, // for Record
          content[0],
          static_cast<uint32_t>(mid),
          static_cast<uint32_t>(end),
          *content[6].begin(),
          "tRF",
          gene_name + (strand == "+" ? "-3p" : "-5p")
        }
      );
    }
    file.close();
  }

  template<class Path>
  bedsix_v read_gencode_trna_gff_into(Path&& input_file_path) {
    bedsix_v results;
    read_gencode_trna_gff(input_file_path, results);
    return results;
  }

  template<class Path>
  void read_GtRNAdb_trna_bed(Path&& input_file_path, bedsix_v& container) {
    /// bed example from GtRNAdb
    /**
     * chr1	4913786	4913856	tRX-Arg-NNN-482-1	1000	+	4913786	4913856	0	1	70,	0,
     * chr1	9553154	9553222	tRX-Und-NNN-49-1	1000	-	9553154	9553222	0	1	68,	0,
     * chr1	34434811	34434883	tRNA-Glu-TTC-2-1	1000	-	34434811	34434883	0	1	72,	0,
     */
    if (!std::filesystem::exists(input_file_path)) {
      std::string file_name = input_file_path;
      throw std::runtime_error("read_GtRNAdb_trna_bed file " + file_name + " does not exist!!");
    }
    std::ifstream file(input_file_path);
    std::string line;
    std::string tips_msg = R"(
Start reading the GtRNAdb bed annotation file.
If any `ERROR` or `Segmentation fault` occurred, 
the reason might be the breaking change of the file format of the annotation.
Please check the URL of GtRNAdb: https://gtrnadb.ucsc.edu/faq.html
)";

    // std::cout << tips_msg << std::endl;
    while(std::getline(file, line)) {
      if (line.empty()) continue;
      if (line.at(0) == '#') continue;
      std::vector<std::string> content = split(line, "\t");
      if (content[0].starts_with("chrUn")) continue;
      std::vector<std::string> split_gene_name_parts = split(content[3], "-");
      if (split_gene_name_parts.size() < 3) {
        throw std::runtime_error("read_GtRNAdb_trna_bed tRNA name " + content[3] + " is not valid!");
      }
      std::string gene_name = split_gene_name_parts[1] + '-' + split_gene_name_parts[2];
      for(int i = 3; i < split_gene_name_parts.size(); ++i){
        gene_name += ("_" + split_gene_name_parts[i]);
      }

      std::string strand = content[5];
      std::size_t start = std::stoul(content[1]); // 0-based (bed)
      std::size_t end = std::stoul(content[2]); // 0-based (bed)
      std::size_t mid = (end - start) / 2 + start;
      container.emplace_back( 
        bedsix {
          {}, // for Record
          content[0],
          static_cast<uint32_t>(start),
          static_cast<uint32_t>(mid),
          strand.front(),
          "tRF",
          gene_name + (strand == "+" ? "-5p" : "-3p")
        }
      );
      container.emplace_back( 
        bedsix {
          {}, // for Record
          content[0],
          static_cast<uint32_t>(mid),
          static_cast<uint32_t>(end),
          strand.front(),
          "tRF",
          gene_name + (strand == "+" ? "-3p" : "-5p")
        }
      );
    }
    file.close();
  }

  template<class Path>
  bedsix_v read_GtRNAdb_trna_bed_into(Path&& input_file_path) {
    bedsix_v results;
    read_GtRNAdb_trna_bed(input_file_path, results);
    return results;
  }
} // namespace bedsixreader

}  // namespace biovoltron
