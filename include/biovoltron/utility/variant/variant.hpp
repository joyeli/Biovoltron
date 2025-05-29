#pragma once

#include <biovoltron/file_io/vcf.hpp>
#include <biovoltron/utility/genotype/genotype.hpp>
#include <biovoltron/utility/range/range_utils.hpp>

namespace biovoltron {

struct Variant {
  Interval location;
  std::string ref;
  std::string alt;
  std::vector<std::string> alleles;
  Genotype gt;
  std::vector<int> pls;
  int gq;
  double qual{};

  auto
  operator==(const Variant& other) const noexcept {
    return std::tie(location, ref, alt)
           == std::tie(other.location, other.ref, other.alt);
  }

  auto
  operator<=>(const Variant& other) const noexcept {
    return std::tie(location, ref, alt)
           <=> std::tie(other.location, other.ref, other.alt);
  }

  auto
  size() const noexcept {
    return location.size();
  }

  auto
  is_snp() const noexcept {
    return ref.size() == alt.size();
  }

  auto
  is_insertion() const noexcept {
    return ref.size() < alt.size();
  }

  auto
  is_deletion() const noexcept {
    return ref.size() > alt.size();
  }

  auto
  to_string() const {
    auto oss = std::ostringstream{};
    oss << location.chrom << "\t" << location.begin + 1 << "\t.\t"
        << (alleles.empty() ? "." : alleles.front()) << "\t";
    auto alts = alleles | std::views::drop(1);
    if (alts.empty())
      oss << ".";
    else
      RangeUtils::format_print(alts, oss);
    oss << "\t" << qual << "\t.\t.\t"
        << "GT:GQ:PL\t" << gt << ":" << gq << ":";
    if (pls.empty())
      oss << ".";
    else
      RangeUtils::format_print(pls, oss);
    return oss.str();
  }

  operator auto() const {
    auto record = VcfRecord{};
    std::istringstream{to_string()} >> record;
    return record;
  }
};

}  // namespace biovoltron