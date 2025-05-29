#pragma once

#include <biovoltron/file_io/core/record.hpp>
#include <biovoltron/utility/interval.hpp>

namespace biovoltron {

/**
 * @ingroup file_io
 * @brief
 * Generic Feature Format is a file format that used to describe genes 
 * or other features of DNA, RNA and protein sequences. 
 * 
 * GFFRecord is a data structure to hold each sample 
 * in the body information of GFF file. 
 * Undefined fields are replaced with "."(dot).
 * 
 * Example:
 * ```cpp
 * #include <iostream>
 * #include <sstream>
 * #include <cassert>
 * #include <biovoltron/utility/interval.hpp>
 * #include <biovoltron/file_io/gff.hpp>
 *
 * int main() {
 *   using namespace biovoltron;
 *
 *   {
 *     GffRecord r;
 *     auto iss = std::stringstream{"ctg123\t.\tmRNA\t10000\t15000\t0\t+\t0\tID=mrna0002;Parent=operon001;Name=subsonicsquirrel"};
 *     iss >> r;
 *     std::cout << r << std::endl;
 *   }
 *
 *   {
 *     GffRecord r1;
 *     GffRecord r2;
 *     auto iss1 = std::stringstream{"ctg123\t.\tmRNA\t10000\t15000\t0\t+\t0\tID=mrna0002;Parent=operon001;Name=subsonicsquirrel"};
 *     iss1 >> r1;
 *     auto iss2 = std::stringstream{"btg123\t.\tmRNA\t10000\t15000\t.\t+\t.\tID=mrna0002;Parent=operon001;Name=subsonicsquirrel"};
 *     iss2 >> r2;
 *
 *     assert(true==(r1<r2));
 *     assert(false==(r1==r2));
 *     assert(false==(r1>r2));
 *   }
 *
 *   {
 *     GffRecord r;
 *     auto iss = std::stringstream{"btg123\t.\tmRNA\t10000\t15000\t.\t+\t.\tID=mrna0002;Parent=operon001;Name=subsonicsquirrel"};
 *     iss >> r;
 *     Interval i = r;
 *   }
 * }
 * ```
 */
struct GffRecord : Record {
  /**
   * @brief The name of sequence where feature is.
   */
  std::string seqid = ".";

  /**
   * @brief Describe source of this feature.
   *
   * Typically the name of software, such as "Genescan" or 
   * name of database, such as "Genebank".
   */
  std::string source = ".";

  /**
   * @brief The type of feature(e.g. gene, DNA, mRNA, exon, and so on).
   *
   * In a well structured GFF file
   * , all the children feature always follow their parent with a single block.
   */
  std::string type = ".";

  /**
   * @brief Start of position of the feature.
   */
  std::uint32_t start{0};

  /**
   * @brief End of position of the feature.
   */
  std::uint32_t end{0};

  /**
   * @brief Score of the feature.
   *
   * A floating number that indicates the confidence of the source in annotated feature,
   */
  float score = 0;

  /**
   * @brief The strand of the feature. 
   *
   * The values of "+" is positive(5'->3'), "-" is negative(3'->5'), and "." is undetermined.
   */
  char strand{'.'};

  /**
   * @brief Phase of CDS features, where CDS means "CoDing Sequence".
   *
   * The phase is one of the integers 0, 1, or 2, indicating the number of bases 
   * that should be removed from the beginning of this feature to reach the first base of the next codon.
   */
  int phase{0};

  /**
   * @brief Providing a list of additional attributes of the feature in the format tag=value. 
   *
   * Multiple pairs are separated by semicolons. 
   * i.e "ID=mrna0002;Parent=operon001;Name=subsonicsquirrel"
   */
  std::string attrs = ".";

  /**
   * @brief Compare two samples in genome feature format.
   *
   * Compare their chromosome sequence first, and then their position.
   */
  auto
  operator<=>(const GffRecord& other) const noexcept {
    return std::tie(seqid, start, end)
           <=> std::tie(other.seqid, other.start, other.end);
  }

  /**
   * @brief Implicit convert to Intervel.
   *
   * Interval is a segment of sequence with name of chromosome, 
   * start and end of position, and the strand of segment.
   */
  operator auto() const { return Interval{seqid, start - 1, end, strand}; }
};

}  // namespace biovoltron
