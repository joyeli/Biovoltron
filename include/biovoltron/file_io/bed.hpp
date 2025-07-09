#pragma once

#include <biovoltron/file_io/core/header.hpp>
#include <biovoltron/file_io/core/record.hpp>
#include <biovoltron/utility/interval.hpp>

namespace biovoltron {

/**
 * This is a data structure header information for genomic features and annotations in the browser extensible data (BED) format.
 * - START_SYMBOLS define as "browser", "track", "#".
 *   - For example: "browser position chr19:49302001-49304701\n" with START_SYMBOLS "browser", we can define this string as header information.
 */
struct BedHeader : Header {
  constexpr static auto START_SYMBOLS = std::array{"browser", "track", "#"};
};

/**
 * This is a data structure for genomic features and annotations in the browser extensible data (BED) format.
 * - As described on the UCSC Genome Browser website.
 */
struct BedRecord : HeaderableRecord {
  BedHeader* header = nullptr;

  /**
   * The name of the chromosome on which the genome feature exists.
   * - Any string can be used. For example, “chr1”, “III”, “myChrom”, “contig1112.23”.
   * - This column is required.
   */
  std::string chrom;

  /**
   * The zero-based starting position of the feature in the chromosome.
   * - The first base in a chromosome is numbered 0.
   * - The start position in each BED feature is therefore interpreted to be 1 greater than the start position listed in the feature.
   *   - For example, start=9, end=20 is interpreted to span bases 10 through 20,inclusive.
   * - This column is required.
   */
  std::uint32_t start{};

  /**
   * The one-based ending position of the feature in the chromosome.
   * - The end position in each BED feature is one-based.
   * - This column is required.
   */
  std::uint32_t end{};

  /**
   * Defines the name of the BED feature.
   * - Any string can be used.
   *   - For example, “LINE”, “Exon3”, “HWIEAS_0001:3:1:0:266#0/1”, or “my_Feature”.
   * - This column is optional.
   */
  std::string name;
  
  /**
   * Defines the score of the BED feature.
   * - BED score range from 0 to 1000
   * - This column is optional.
   */
  std::int32_t score = 0;

  /**
   * Defines the the strand.
   * - "+" : sense
   * - "-" : anti-sense
   * - This column is optional.
   */
  char strand = 0;

  /**
   * The starting position at which the feature is drawn thickly.
   * - Coding sequence(CDS) start position
   * - Allowed yet ignored by bedtools.
   */
  std::uint32_t thick_start = 0;

  /**
   * The ending position at which the feature is drawn thickly.
   * - Coding sequence(CDS) end position
   * - Allowed yet ignored by bedtools.
   */
  std::uint32_t thick_end = 0;

  /**
   * An RGB value of the form R,G,B.
   * - e.g. 255,0,0
   * - Allowed yet ignored by bedtools.
   */
  std::string item_rgb = "0,0,0";

  /**
   * The number of blocks (exons) in the BED line.
   * - pair with block_sizes and block_starts 
   * - Allowed yet ignored by bedtools.
   */
  std::uint32_t block_count = 0;

  /**
   * A comma-separated list of the block sizes.
   * - pair with block_count
   *   - If block_count = 3 then block_sizes will be 3 number composite by ",".
   *   - For example: block_sizes = "354,109,1189"
   * - Allowed yet ignored by bedtools.
   */
  std::string block_sizes = "0";

  /**
   * A comma-separated list of block starts.
   * - pair with block_count
   *   - If block_count = 3 then block_starts will be 3 number composite by ",".
   *   - For example: block_starts = "0,739,1347"
   * - Allowed yet ignored by bedtools.
   */
  std::string block_starts = "0";

  /**
   * Three-way comparation (spaceship operator): strong ordering
   * - compare two BedRecord order with chrom->start->end
   *   - less:
   *     - if this{“chr1”, 31, 500}, other{“chr1”, 36, 560}
   *     - then this <=> other < 0 because “chr1” == “chr1” and 31 < 36
   *   - equal:
   *     - if this{“chr1”, 31, 500}, other{“chr1”, 31, 500}
   *     - then this <=> other == 0 because “chr1” == “chr1”, 31 == 31 and 500 == 500
   *   - greater:
   *     - if this{“chr1”, 39, 500}, other{“chr1”, 36, 560}
   *     - then this <=> other > 0 because “chr1” == “chr1” and 39 > 36
   */
  auto
  operator<=>(const BedRecord& other) const noexcept {
    return std::tie(chrom, start, end)
       <=> std::tie(other.chrom, other.start, other.end);
  }

  /**
   * Implicitly converted to Interval
   * - return a Interval include chrom, start, end, strand.
   */
  operator auto() const { return Interval{chrom, start, end, strand}; }
};

/**
 * This is a data structure for drawing BedGraph Track Format as graph.
 * - The bedGraph format allows display of continuous-valued data in track format.
 * - This display type is useful for probability scores and transcriptome data.
 * - Data structure combine BED and WIG format:
 *   - only part of BED and WIG format
 */
struct BedGraphRecord : HeaderableRecord {
  BedHeader* header = nullptr;

  /**
   * The name of the chromosome on which the genome feature exists.
   * - Any string can be used. For example, “chr1”, “III”, “myChrom”, “contig1112.23”.
   */
  std::string chrom;

  /**
   * The starting position of the feature in the chromosome.
   */  
  int start{};

  /**
   * The ending position of the feature in the chromosome.
   */  
  int end{};

  /**
   * Defines the score of the annotations.
   * - score can be GC percent, probability scores, and transcriptome data.
   * - BedGraph track data values can be integer or real, positive or negative values.
   */  
  float score{};

  /**
   * Three-way comparation (spaceship operator): strong ordering
   * - compare two BedRecord order with chrom->start->end
   *   - less:
   *     - if this{“chr1”, 31, 500}, other{“chr1”, 36, 560}
   *     - then this <=> other < 0 because “chr1” == “chr1” and 31 < 36
   *   - equal:
   *     - if this{“chr1”, 31, 500}, other{“chr1”, 31, 500}
   *     - then this <=> other == 0 because “chr1” == “chr1”, 31 == 31 and 500 == 500
   *   - greater:
   *     - if this{“chr1”, 39, 500}, other{“chr1”, 36, 560}
   *     - then this <=> other > 0 because “chr1” == “chr1” and 39 > 36
   */
  auto
  operator<=>(const BedRecord& other) const noexcept {
    return std::tie(chrom, start, end)
       <=> std::tie(other.chrom, other.start, other.end);
  }
};

}  // namespace biovoltron
