#pragma once

#include <limits>
#include <stdexcept>
#include <string>

namespace biovoltron {

struct Interval {
  static constexpr auto CHROM_SEPARATOR = ':';
  static constexpr auto BEGIN_END_SEPARATOR = '-';
  static constexpr auto END_OF_CHROM = '+';
  static constexpr auto DIGIT_SEPARATOR = ',';

  std::string chrom;
  std::uint32_t begin{};
  std::uint32_t end{};
  char strand = '+';

 private:
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
  Interval() = default;
  Interval(std::string chrom, std::uint32_t begin, std::uint32_t end, char strand = '+')
  : chrom(std::move(chrom)), begin(begin), end(end), strand(strand){
    if (strand != '+' && strand != '-')
      throw std::invalid_argument("invalid strand symbol");
    if (end < begin)
      throw std::invalid_argument("invalid end must not be less than begin");
  }

  Interval(
    std::convertible_to<const std::string&> auto const& interval_string) {
    parse_string(interval_string);
    if (end < begin)
      throw std::invalid_argument("invalid interval string");
  }

  auto
  size() const noexcept {
    return end - begin;
  }

  auto
  empty() const noexcept {
    return size() == 0;
  }

  auto
  overlaps(const Interval& other) const noexcept {
    return chrom == other.chrom && strand == other.strand && begin < other.end && other.begin < end;
  }

  auto
  contains(const Interval& other) const noexcept {
    return chrom == other.chrom && strand == other.strand && begin <= other.begin && end >= other.end;
  }

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

  auto
  expand_with(std::uint32_t padding) const {
    return Interval{chrom, begin - padding, end + padding, strand};
  }

  auto
  to_string() const {
    return strand + chrom + CHROM_SEPARATOR + std::to_string(begin) + BEGIN_END_SEPARATOR
         + std::to_string(end);
  }

  auto
  operator<=>(const Interval&) const noexcept = default;
};

}  // namespace biovoltron
