#pragma once

#include <limits>
#include <stdexcept>
#include <string>

namespace biovoltron {

/**
 * @ingroup utility
 * @brief Genomic interval with chromosome, range, and strand.
 *
 * Example:
 * ```cpp
 * #include <iostream>
 * #include <biovoltron/utility/interval.hpp>
 * 
 * int main() {
 *     using namespace biovoltron;
 *     Interval iv1{"chr1", 100, 200};       // + strand by default
 *     Interval iv2{"-chr2:500-600"};        // parse from string
 *     if (iv1.overlaps(iv2)) {
 *         std::cout << "Intervals overlap.\n";
 *     }
 *     else {
 *         std::cout << "Intervals do not overlap.\n";
 *     }
 * }
 * ```
 */
struct Interval {

  /**
   * Separators for genomic intervals.
   * - CHROM_SEPARATOR: Separator between chromosome and positions (e.g., "chr1:100-200").
   * - BEGIN_END_SEPARATOR: Separator between begin and end positions.
   * - END_OF_CHROM: Symbol for the end of chromosome.
   * - DIGIT_SEPARATOR: Digit separator (ignored during parsing).
   */
  static constexpr auto CHROM_SEPARATOR = ':';
  static constexpr auto BEGIN_END_SEPARATOR = '-';
  static constexpr auto END_OF_CHROM = '+';
  static constexpr auto DIGIT_SEPARATOR = ',';

  /**
   * Chromosome name.
   */
  std::string chrom;
  /**
   * Start position of the interval (inclusive, 0-based).
   */
  std::uint32_t begin{};
  /**
   * End position of the interval (exclusive, 0-based).
   */
  std::uint32_t end{};
  /**
   * Strand of the interval, either '+' or '-', Default is '+'.
   */
  char strand = '+';

 private:
  /**
   * Parse a string into an interval.
   *
   * @param interval_string Interval in string form.
   */
  auto
  parse_string(const std::string& interval_string) {
    auto iv_str = interval_string;
    if (interval_string.front() == '+' || interval_string.front() == '-') {
      strand = interval_string.front();
      iv_str = iv_str.substr(1);
    } 

    if (const auto colon = iv_str.find(CHROM_SEPARATOR);
        colon == std::string::npos) {
      chrom = iv_str;
      begin = 0;
      end = std::numeric_limits<std::uint32_t>::max();
    } else {
      chrom = iv_str.substr(0, colon);
      auto remain = iv_str.substr(colon + 1);
      std::erase(remain, DIGIT_SEPARATOR);
      begin = std::stoul(remain);
      if (const auto dash = remain.find(BEGIN_END_SEPARATOR);
          dash == std::string::npos) {
        if (remain.back() == END_OF_CHROM)
          end = std::numeric_limits<std::uint32_t>::max();
        else
          end = begin + 1;
      } else
        end = std::stoul(remain.substr(dash + 1));
    }
  }

 public:
  /**
   * Default constructor.
   */
  Interval() = default;
  /**
   * Constructor with chromosome, begin, end, and strand.
   *
   * @param chrom Chromosome name.
   * @param begin Start position (inclusive).
   * @param end End position (exclusive).
   * @param strand Strand of the interval ('+' or '-').
   * @throw std::invalid_argument If strand is invalid or end < begin.
   */
  Interval(std::string chrom, std::uint32_t begin, std::uint32_t end, char strand = '+')
  : chrom(std::move(chrom)), begin(begin), end(end), strand(strand){
    if (strand != '+' && strand != '-')
      throw std::invalid_argument("invalid strand symbol");
    if (end < begin)
      throw std::invalid_argument("invalid end must not be less than begin");
  }

  /**
   * Constructor from a string representation of an interval.
   *
   * @param interval_string Interval in string form (e.g., "+chr1:100-200").
   * @throw std::invalid_argument If the interval string is invalid.
   */
  Interval(
    std::convertible_to<const std::string&> auto const& interval_string) {
    parse_string(interval_string);
    if (end < begin)
      throw std::invalid_argument("invalid interval string");
  }

  /**
   * Get interval length.
   * 
   * @return Length of the interval (end - begin).
   */
  auto
  size() const noexcept {
    return end - begin;
  }

  /**
   * Checks if the interval is empty.
   * 
   * @return true if the interval is empty, false otherwise
   */
  auto
  empty() const noexcept {
    return size() == 0;
  }

  /**
   * Checks if the interval overlaps with another interval.
   * 
   * @param other The other interval to check against.
   * @return true if the intervals overlap, false otherwise
   */
  auto
  overlaps(const Interval& other) const noexcept {
    return chrom == other.chrom && strand == other.strand && begin < other.end && other.begin < end;
  }

  /**
   * Checks if the interval contains another interval.
   * 
   * @param other The other interval to check against.
   * @return true if this interval contains the other interval, false otherwise
   */
  auto
  contains(const Interval& other) const noexcept {
    return chrom == other.chrom && strand == other.strand && begin <= other.begin && end >= other.end;
  }

  /**
   * Return a new interval that spans both this interval and another interval.
   * @param other Other interval.
   * @return A new interval spanning both.
   * @throw std::invalid_argument If chromosomes or strands differ.
   */
  auto
  span_with(const Interval& other) const {
    if (chrom != other.chrom)
      throw std::invalid_argument(
        "Interval::span_with(): Cannot get span for intervals on different "
        "chroms.");
    if (strand != other.strand)
      throw std::invalid_argument(
        "Interval::span_with(): Cannot get span for intervals on different "
        "strands.");
    return Interval{chrom, std::min(begin, other.begin),
                    std::max(end, other.end), strand};
  }

  /**
   * Expand the interval by a specified padding.
   * 
   * @param padding The amount to expand the interval on both sides.
   * @return A new interval expanded by the specified padding.
   */
  auto
  expand_with(std::uint32_t padding) const {
    return Interval{chrom, begin - padding, end + padding, strand};
  }

  /**
   * Convert the interval to a string representation.
   * 
   * @return A string representation of the interval.
   */
  auto
  to_string() const {
    return strand + chrom + CHROM_SEPARATOR + std::to_string(begin) + BEGIN_END_SEPARATOR
         + std::to_string(end);
  }

  /**
   * Comparison operator for intervals.
   */
  auto
  operator<=>(const Interval&) const noexcept = default;
};

}  // namespace biovoltron
