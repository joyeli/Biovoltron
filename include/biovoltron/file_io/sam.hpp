#pragma once

#include <biovoltron/file_io/cigar.hpp>
#include <biovoltron/file_io/core/header.hpp>
#include <biovoltron/file_io/core/record.hpp>
#include <biovoltron/utility/interval.hpp>
#include <biovoltron/utility/read/quality_utils.hpp>

namespace biovoltron {

/**
 * @ingroup file_io
 * @brief A class about the config of SAM file and some core methods
 */
struct SamUtil {
  /**
   * @brief Maximum length of reads.
   */
  static constexpr auto MAX_READ_LENGTH = 256;

  /**
   * @brief Gap open penalty defined as 40. Represent as quality value characters as in fastq file.
   */
  static inline const auto GAP_OPEN_PENALTY =
    std::string(MAX_READ_LENGTH, 40 + QualityUtils::ASCII_OFFSET);

  /**
   * @brief Gap continuation penalty defined as 10. Represent as quality value characters as in fastq file
   */
  static inline const auto GAP_CONTINUATION_PENALTY =
    std::string(MAX_READ_LENGTH, 10 + QualityUtils::ASCII_OFFSET);

  /**
   * @brief Bitwise flags of a read alignment.
   * - READ_PAIRED: read is paired
   * - PROPER_PAIR: read mapped in proper pair
   * - READ_UNMAPPED: read unmapped
   * - MATE_UNMAPPED: mate unmapped
   * - READ_REVERSE_STRAND: read is on the reverse strand
   * - MATE_REVERSE_STRAND: mate is on the reverse strand
   * - FIRST_OF_PAIR: this is read1
   * - SECOND_OF_PAIR: this is read2
   * - SECONDARY_ALIGNMENT: this is not a primary alignment
   * - READ_FAILS_QUALITY_CHECK: alignment failed the quality check
   * - DUPLICATE_READ: this is a PCR or optical duplicate
   * - SUPPLEMENTARY_ALIGNMENT: this is a supplementary alignment
   */
  enum Flag {
    READ_PAIRED = 0x1,
    PROPER_PAIR = 0x2,
    READ_UNMAPPED = 0x4,
    MATE_UNMAPPED = 0x8,
    READ_REVERSE_STRAND = 0x10,
    MATE_REVERSE_STRAND = 0x20,
    FIRST_OF_PAIR = 0x40,
    SECOND_OF_PAIR = 0x80,
    SECONDARY_ALIGNMENT = 0x100,
    READ_FAILS_QUALITY_CHECK = 0x200,
    DUPLICATE_READ = 0x400,
    SUPPLEMENTARY_ALIGNMENT = 0x800
  };

  /**
   * @brief Orientation of the read.
   * - FR: forward-reverse
   * - FF: forward-forward
   * - RR: reverse-reverse
   * - RF: reverse-forward
   */
  enum Orientation { FR, FF, RR, RF };

  /**
   * @brief Returns the orientation of the read.
   * 
   * @param read_forward whether the read is mapped to forward strand
   * @param mate_forward whether the mate is mapped to forward strand
   * @return The orientation of the read.
   */
  static auto
  compute_ori(bool read_forward, bool mate_forward) noexcept {
    if (read_forward != mate_forward)
      return read_forward ? FR : RF;
    return read_forward ? FF : RR;
  }

  /**
   * @brief Returns the length of the template.
   * 
   * @param read_pos position of the read
   * @param read_cigar CIGAR of the read
   * @param read_forward whether the read is mapped to forward strand
   * @param mate_pos position of the mate
   * @param mate_cigar CIGAR of the mate
   * @param mate_forward whether the mate is mapped to forward strand
   * @return The length which the fragment reference to.
   */
  static auto
  compute_tlen(std::int32_t read_pos, const Cigar& read_cigar,
               bool read_forward, std::int32_t mate_pos,
               const Cigar& mate_cigar, bool mate_forward) noexcept
    -> std::int32_t {
    if (read_pos > mate_pos)
      return -compute_tlen(mate_pos, mate_cigar, mate_forward, read_pos,
                           read_cigar, read_forward);

    switch (const auto read_ori = compute_ori(read_forward, mate_forward);
            read_ori) {
      case FR:
        return mate_pos + mate_cigar.ref_size() - read_pos;
      case FF:
        if (const auto tlen = mate_pos + mate_cigar.read_size()
                            - (read_pos + read_cigar.read_size());
            tlen != 0)
          return tlen + (tlen > 0 ? 1 : -1);
        else
          return 0;
      case RR:
        if (const auto tlen = mate_pos + mate_cigar.ref_size()
                            - (read_pos + read_cigar.ref_size());
            tlen != 0)
          return tlen + (tlen > 0 ? 1 : -1);
        else
          return 0;
      case RF:
        if (const auto tlen = mate_pos - (read_pos + read_cigar.ref_size()) + 1;
            tlen != 0)
          return tlen + (tlen > 0 ? 1 : -1);
        else
          return 0;
      default:
        return 0;
    }
  }
};

/**
 * @ingroup file_io
 * @brief A class about the header of SAM format
 */
struct SamHeader : Header {
  /**
   * @brief The indicator of sam file header
   */
  constexpr static auto START_SYMBOLS = std::array{"@"};
};

/**
 * @tparam Encoded
 * 
 * @ingroup file_io
 * @brief A class that save the alignment section of SAM file
 * 
 * Example
 * ```cpp
 * #include <iostream>
 * #include <biovoltron/file_io/sam.hpp>
 * 
 * int main()
 * {
 *     using namespace biovoltron;
 *     auto record1 =
 *         SamRecord{
 *             {},
 *             nullptr,
 *             "r1",
 *             SamUtil::READ_PAIRED + SamUtil::PROPER_PAIR + SamUtil::MATE_REVERSE_STRAND + SamUtil::FIRST_OF_PAIR,
 *             "ref",
 *             2,
 *             30,
 *             "3M1D2M2I6M",
 *             "=",
 *             22,
 *             29,
 *             "TCGACGGTGACTG",
 *             "*",
 *             {}};
 * 
 *     if (record1.read_paired())
 *         std::cout << "record1 is paired.\n";
 *     if (record1.proper_pair())
 *         std::cout << "record1 is in proper pair.\n";
 *     std::cout << "record1 starts from position " << record1.pos << " with a size of " << record1.size() << ".\n";
 * }
 * ```
 */
template<bool Encoded = false>
struct SamRecord : HeaderableRecord {
  /**
   * @brief Whether this file is encoded.
   */
  constexpr static auto encoded = Encoded;

  /**
   * @brief Pointer to the header.
   */
  SamHeader* header = nullptr;

  /**
   * @brief Name of the aligned read.
   */
  std::string qname;

  /**
   * @brief Bitwise flag of the aligned read.
   */
  std::uint16_t flag{};

  /**
   * @brief Reference sequence name of the alignment.
   */
  std::string rname;

  /**
   * @brief The leftmost mapping position of the first matching base.
   */
  std::uint32_t pos{};

  /**
   * @brief Mapping quality.
   */
  std::uint16_t mapq{};

  /**
   * @brief CIGAR describe the alignment between read and reference.
   */
  Cigar cigar;

  /**
   * @brief Reference sequence name of the mate read. 
   * Set to "=" if equal to rname, "*" if the information is unavailable.
   */
  std::string rnext;

  /**
   * @brief The leftmost mapping postion of the first matching base of the mate read. Set to 0 if is single-end read.
   */
  std::uint32_t pnext{};

  /**
   * @brief The length of the template. 
   * Positive sign for the leftmost segment and negative sign for the rightmost segment.The sign of the segments in the middle is undefined.
   * Set to 0 for single-end read or when the information is unavailable.
   */
  std::int32_t tlen{};

  /**
   * @brief Sequence of the read. Can be "*" if the sequence is not stored.
   */
  std::conditional_t<Encoded, istring, std::string> seq;

  /**
   * @brief ASCII of quality score. Can be "*" if quality is not stored.
   */
  std::string qual;

  /**
   * @brief Option field. Stored as TAG:TYPE:VAL.
   */
  std::vector<std::string> optionals;

  /**
   * @brief Checks if the read is paired.
   * 
   * @return true if the read is paired, false otherwise
   */
  auto
  read_paired() const noexcept {
    return !!(flag & SamUtil::READ_PAIRED);
  }

  /**
   * @brief Checks if the read is mapped in proper pair.
   * 
   * @return true if the read is mapped in proper pair, false otherwise
   */
  auto
  proper_pair() const noexcept {
    return !!(flag & SamUtil::PROPER_PAIR);
  }

  /**
   * @brief Checks if the read is unmapped.
   * 
   * @return true if the read is unmapped, false otherwise
   */
  auto
  read_unmapped() const noexcept {
    return !!(flag & SamUtil::READ_UNMAPPED);
  }

  /**
   * @brief Checks if the mate is unmapped.
   * 
   * @return true if the mate is unmapped, false otherwise
   */
  auto
  mate_unmapped() const noexcept {
    return !!(flag & SamUtil::MATE_UNMAPPED);
  }

  /**
   * @brief Checks if the read is on the reverse strand.
   * 
   * @return true if the read is on the reverse strand, false otherwise
   */
  auto
  read_reverse_strand() const noexcept {
    return !!(flag & SamUtil::READ_REVERSE_STRAND);
  }

  /**
   * @brief Checks if the mate is on the reverse strand.
   * 
   * @return true if the mate is on the reverse strand, false otherwise
   */
  auto
  mate_reverse_strand() const noexcept {
    return !!(flag & SamUtil::MATE_REVERSE_STRAND);
  }

  /**
   * @brief Checks if the read is the first read.
   * 
   * @return true if the read is the first read, false otherwise
   */
  auto
  first_of_pair() const noexcept {
    return !!(flag & SamUtil::FIRST_OF_PAIR);
  }

  /**
   * @brief Checks if the read is the second read.
   * 
   * @return true if the read is the second read, false otherwise
   */
  auto
  second_of_pair() const noexcept {
    return !!(flag & SamUtil::SECOND_OF_PAIR);
  }

  /**
   * @brief Checks if the alignment is a secondary alignment.
   * 
   * @return true if the alignment is a secondary alignment, false otherwise
   */
  auto
  secondary_alignment() const noexcept {
    return !!(flag & SamUtil::SECONDARY_ALIGNMENT);
  }

  /**
   * @brief Checks if the alignment failed the quality check.
   * 
   * @return true if the alignment failed the qualtiy check, false otherwise
   */
  auto
  read_fails_quality_check() const noexcept {
    return !!(flag & SamUtil::READ_FAILS_QUALITY_CHECK);
  }

  /**
   * @brief Checks if the read is a PCR or optical duplicate.
   * 
   * @return true if the read is a PCR or optical duplicate, false otherwise
   */
  auto
  duplicate_read() const noexcept {
    return !!(flag & SamUtil::DUPLICATE_READ);
  }

  /**
   * @brief Checks if the alignment is a supplementary alignment.
   * 
   * @return true if the alignment is a supplementary alignment, false otherwise
   */  
  auto
  supplementary_alignment() const noexcept {
    return !!(flag & SamUtil::SUPPLEMENTARY_ALIGNMENT);
  }

  /**
   * @brief Returns the number of the bases in the read.
   * 
   * @return the number of the bases in the read
   */
  auto
  size() const noexcept {
    return seq.size();
  }

  /**
   * @brief Checks if the sequence of the read is empty.
   * 
   * @return true if the sequence of the read is empty, false otherwise
   */
  auto
  empty() const noexcept {
    return seq.empty();
  }

  /**
   * @brief Returns the 0-base position of the first matching base(i.e. pos - 1).
   * 
   * @return the 0-bae position for the first matching base
   */
  auto
  begin() const noexcept {
    return pos - 1;
  }

  /**
   * @brief Returns the 0-base position of the last matching base.
   * 
   * @return the 0-base postion of the last matching base
   */
  auto
  end() const noexcept {
    return begin() + cigar.ref_size();
  }

  /**
   * @brief 
   * 
   * @return tlen - 1  
   */
  auto
  mate_begin() const noexcept {
    return tlen - 1;
  }

  auto
  tlen_well_defined() const noexcept {
    if (tlen == 0)
      return false;
    if (!read_paired())
      return false;
    if (read_unmapped() || mate_unmapped())
      return false;
    if (read_reverse_strand() == mate_reverse_strand())
      return false;
    if (read_reverse_strand())
      return end() > mate_begin() + 1;
    return begin() <= mate_begin() + tlen;
  }

  /**
   * @brief Generate a string with length same as seq and fill with GAP_OPEN_PENALTY.
   * 
   * @return a string with length same as seq and fill with GAP_OPEN_PENALTY
   */
  auto
  insertion_gop() const noexcept {
    return std::string_view{SamUtil::GAP_OPEN_PENALTY.data(), seq.size()};
  }

  /**
   * @brief Generate a string with length same as seq and fill with GAP_OPEN_PENALTY.
   * 
   * @return a string with length same as seq and fill with GAP_OPEN_PENALTY
   */
  auto
  deletion_gop() const noexcept {
    return std::string_view{SamUtil::GAP_OPEN_PENALTY.data(), seq.size()};
  }

  /**
   * @brief Generate a string with length same as seq and fill with GAP_CONTINUATION_PENALTY.
   * 
   * @return a string with length same as seq and fill with GAP_CONTIUATION_PENALTY
   */
  auto
  overall_gcp() const noexcept {
    return std::string_view{SamUtil::GAP_CONTINUATION_PENALTY.data(),
                            seq.size()};
  }

  /**
   * @brief Compares with other. First compares rname to other.rname, then pos ot other.pos
   * 
   * @param other SamRecord to compare with
   * @return true if the corresponding comparison holds, false otherwise.
   */
  auto
  operator<=>(const SamRecord& other) const noexcept {
    return std::tie(rname, pos) <=> std::tie(other.rname, other.pos);
  }

  /**
   * @brief Returns a struct Interval of {chrom = rname, begin = begin(), end = end(), strand = '-'} if read_reverse_strand() == true,
   * else returns a struct Inverval of {chrom = rname, begin = begin(), end = end(), strand = '+'}
   * 
   * @return Interval{rname, begin(), end(), '-'} if read_reverse_strand() == true,
   *         Interval{rname, begin(), end(), '+'} otherwise
   */
  operator auto() const { 
    auto strand = read_reverse_strand() ? '-' : '+';
    return Interval{rname, begin(), end(), strand}; 
  }
};

}  // namespace biovoltron
