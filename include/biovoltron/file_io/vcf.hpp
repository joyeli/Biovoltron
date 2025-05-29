#pragma once

#include <biovoltron/file_io/core/header.hpp>
#include <biovoltron/file_io/core/record.hpp>
#include <biovoltron/utility/interval.hpp>

namespace biovoltron {
/**
 * @ingroup file_io
 * @brief
 * VCF file: Files of variant calling format store only variantions 
 * along with reference genome. 
 * 
 * VcfHeader is a data structure to store header of VCF file.
 *
 * The header begins the file and provides metadata describing the body of the file. 
 * Lines of header are denoted as starting with #. 
 * 
 * Example:
 * ```cpp
 * #include <iostream>
 * #include <sstream>
 * #include <biovoltron/utility/interval.hpp>
 * #include <biovoltron/file_io/vcf.hpp>
 *
 * int main() {
 *   using namespace biovoltron;
 *
 *   VcfHeader h;
 *   auto iss = std::istringstream{"#This is a line of header.\n"};
 *
 *   ss >> h;
 *   std::cout << h << std::endl;
 * }
 * ```
 */
struct VcfHeader : Header {
  constexpr static auto START_SYMBOLS = std::array{"#"};
};

/**
 * @ingroup file_io
 * @brief
 * VCF file: Files of variant calling format store only variantions 
 * along with reference genome. 
 *
 * A data structure to hold each sample 
 * in the body information of VCF file.
 *
 * For each sample, separate each field with a tab("\t").
 *
 * Example:
 * ```cpp
 * #include <iostream>
 * #include <sstream>
 * #include <cassert>
 * #include <biovoltron/utility/interval.hpp>
 * #include <biovoltron/file_io/vcf.hpp>
 *
 * int main() {
 *   using namespace biovoltron;
 *
 *   {
 *   	VcfRecord r;
 *   	auto iss = std::istringstream{"20\t1110696\trs6040355\tA\tG,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4"};
 *   	iss >> r;
 *   	std::cout << r << std::endl;
 *   }
 *
 *   {
 *   	VcfRecord r1;
 *   	VcfRecord r2;
 *   	auto iss1 = std::stringstream{"20\t1110695\trs6040355\tA\tG,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4"};
 *   	iss1 >> r1;
 *   	auto iss2 = std::stringstream{"20\t1110696\trs6040355\tA\tG,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4"};
 *   	iss2 >> r2;
 *
 *   	assert(true==(r1<r2));
 *   	assert(false==(r1==r2));
 *   	assert(false==(r1>r2));
 *   }
 *
 *   {
 *   	VcfRecord r;
 *   	auto iss = std::stringstream{"20\t1110696\trs6040355\tA\tG,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4"};
 *   	iss >> r;
 *   	Interval i = r;
 *   }
 * }
 * ```
 */
struct VcfRecord : HeaderableRecord {
  /**
   * @brief A pointer to its header with description of corresponding vcf file.
   */
  VcfHeader* header = nullptr;

  /**
   * @brief Name of chromosome on which variation being called.
   */
  std::string chrom{};
  /**
   * @brief Position of variation on given chromosome.
   */
  std::uint32_t pos{};
  /**
   * @brief The identifier of variation, e.g. identifier in dbSNP.
   *
   * If there is no identifier availanble, id is denoted as ".".
   */
  std::string id;
  /**
   * @brief The reference base at given position on reference sequence.
   *
   * For example, the reference allele is "A" at given position of this sample.
   */
  std::string ref;
  /**
   * @brief The possible alternative alleles at smae position of this record.
   *
   * For example, the allele of this record is "C", and the possible alternative alleles are "G,T" at same position.
   */
  std::string alt;
  /**
   * @brief Quality score associated with the sequence of the given alleles.
   * High qualityscores indicate high confidence calls. 
   */
  double qual{};
  /**
   * @brief A flag that indicated a given variation has passed or failed of a set of filters.
   *
   * The filters apply to the sample are described in the header. 
   *
   * i.e. 
   * ```
   * #FILTER=<ID=q10,Description="Quality below 10">
   * #FILTER=<ID=s50,Description="Less than 50% of samples have data">
   * ```
   *
   * Denoted as PASS, if all the filter were passed successfully.
   */
  std::string filter;
  /**
   * @brief An extensible list of key-value pairs of describing the variation.
   *
   * The fields of info are described in the header. 
   *
   * i.e. 
   * ```
   * #INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
   * #INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
   * ```
   *
   * The value at the field would be "NS=3;DP=14" for example.
   */
  std::string info;
  /**
   * @brief An optional extensible list of fields for descibing the samples(reads).
   * 
   * The fields of format are described in the header. 
   *
   * i.e. 
   * ```
   * #FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
   * #FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
   * ```
   */
  std::string format;
  /**
   * @brief Optional: it may be empty.
   *
   * For each sample described in the file, 
   * values are given for the feilds listed in format field.
   *
   * For example, th format filed is "GT:GQ:DP:HQ", 
   * and the value of this field would be "0|0:48:1:51,51".
   *
   * First number is the allele value. The allele values are 0 for the reference allele
   * , 1 for the first allele listed in ALT, 2 for the second allele list in ALT and so on.
   *
   * "|" means genotype phased
   * "/" means genotype unphased
   */
  std::vector<std::string> samples;
  
  /**
   * @brief Compare two samples in variant call format.
   *
   * Compare their chromosome sequence first, and then their position.
   */
  auto
  operator<=>(const VcfRecord& other) const noexcept {
    return std::tie(chrom, pos) <=> std::tie(other.chrom, other.pos);
  }

  /**
   * @brief Implicit convert to Intervel.
   *
   * Interval is a segment of sequence with name of chromosome, 
   * start and end of position, and the strand of segment.
   */
  operator auto() const { return Interval{chrom, pos - 1, pos, '+'}; }
};

}  // namespace biovoltron
