#pragma once

#include <biovoltron/utility/interval.hpp>
#include <biovoltron/utility/istring.hpp>
#include <biovoltron/file_io/sam.hpp>
#include <cassert>
#include <numeric>
#include <set>
#include <vector>
#include <utility>

namespace biovoltron {

struct Mismatch {
  std::uint32_t pos;
  char correct_base;
  bool operator==(const Mismatch&) const = default;
};

struct Hit {
  // The mismatch positions on the original read with the correct base.
  std::vector<Mismatch> mismatches{};
  // The positions on the original read where T to C conversion happened.
  std::set<std::uint32_t> tc_set{};
  Interval intv;
  // ranges is only used when extending toward the 5' end.
  std::vector<std::pair<std::uint32_t, std::uint32_t>> ranges{};
};

struct Alignment {
  std::string name;
  std::string seq;
  std::string qual;
  bool forward;
  std::uint32_t tail_pos{};   // -1 for no tail.
  std::vector<Hit> hits{};
  std::uint32_t counts{};
  // head_pos is only used when extending toward the 5' end.
  std::uint32_t head_pos{};   // -1 for no head.
};

// Note: If map to reverse strand:
// - reverse complement SEQ, TL
// - reverse CIGAR, QUAL, MD
inline auto aln_to_sam_list(const Alignment& aln) {
  if (aln.hits.empty())
    return std::vector<SamRecord<>>{};

  auto tail_len = (aln.tail_pos == -1) ? 0 : aln.seq.size() - aln.tail_pos;
  auto has_mismatch = !aln.hits.front().mismatches.empty();

  auto sam = SamRecord{};
  sam.qname = aln.name;
  sam.flag = aln.forward ? 0 : 16;
  sam.mapq = 255 - tail_len; 
  sam.cigar = (tail_len == 0)
    ? std::to_string(aln.seq.size()) + 'M'
    : std::to_string(aln.tail_pos) + 'M' + std::to_string(tail_len) + 'S';
  if (!aln.forward) std::ranges::reverse(sam.cigar);
  sam.rnext = "*";
  sam.pnext = 0;
  sam.tlen = 0;
  sam.seq = aln.forward ? aln.seq : Codec::rev_comp(aln.seq);
  sam.qual = aln.qual;
  if (!aln.forward) std::ranges::reverse(sam.qual);

  // NH tag: Number of reported alignments.
  sam.optionals.emplace_back("NH:i:" + std::to_string(aln.hits.size()));

  // Tail tag: Report tail.
  if (tail_len != 0) {
    auto tail = aln.forward ? aln.seq.substr(aln.tail_pos) 
      : Codec::rev_comp(aln.seq.substr(aln.tail_pos));
    sam.optionals.emplace_back("TL:Z:" + tail);
  }

  // MD tag: Encoding mismatched and deleted reference bases.
  // Used to reconstruct the mapped reference sequence without
  // requiring access to the reference.
  //
  // [0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*
  //
  // ex: 10A5^AC6 
  // 10 matches followed by a mismatch with 'A' on the reference,
  // and then 5 bases match and 2 base deletion ('AC' is the deltd::vector<Interval> alns{};
  //   std::vector<Mismatch> mismatches{};
  //   ted
  // bases in the reference), and 6 more matches.
  //
  // Softclip (tail) doesn't count.

  // Calculating match segments length between mismatches for MD tag.
  auto match_segs = std::vector<std::uint32_t>{};
  if (has_mismatch) {
    std::ranges::transform(
      aln.hits.front().mismatches,
      std::back_inserter(match_segs),
      [](const auto& mismatch) { return mismatch.pos; }
    );
    match_segs.push_back(aln.seq.size() - tail_len); // Exclude tail.
    std::ranges::sort(match_segs);
    std::adjacent_difference(
      match_segs.begin(), 
      match_segs.end(), 
      match_segs.begin(), 
      [](auto a, auto b) { return a - b - 1; }
    );
  }

  auto sams = std::vector<SamRecord<>>{aln.hits.size(), sam};

  for (auto i = 0; i < aln.hits.size(); i++) {
    sams[i].rname = aln.hits[i].intv.chrom;
    sams[i].pos = aln.hits[i].intv.begin + 1; // 1-based coordinate.
    
    if (has_mismatch) {
      auto mismatches = aln.hits[i].mismatches;
      std::ranges::sort(mismatches, {}, &Mismatch::pos);

      auto md_str = std::string{};
      for (auto i = 0; i < mismatches.size(); i++) {
        md_str += std::to_string(match_segs[i]) 
               +  mismatches[i].correct_base;
      }
      md_str += std::to_string(match_segs.back());

      if (!aln.forward) std::ranges::reverse(md_str);
      sams[i].optionals.emplace_back("MD:Z:" + md_str); 
    }
  }

  return sams;
}

}  // namespace biovoltron
