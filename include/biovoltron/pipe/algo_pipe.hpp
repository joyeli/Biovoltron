#pragma once

#include <biovoltron/algo/align/exact_match/fm_index.hpp>
#include <biovoltron/applications/burrow_wheeler_aligner/burrow_wheeler_aligner.hpp>
#include <biovoltron/applications/haplotypecaller/haplotypecaller.hpp>

namespace biovoltron {

template<int SA_INTV, typename size_type, SASorter Sorter>
auto
operator|(FastaRecord<true> ref, FMIndex<SA_INTV, size_type, Sorter> index) {
  // replace 'N' to 'A'
  std::ranges::replace(ref.seq, 4, 0);
  index.build(ref.seq);
  return index;
}

template<std::ranges::range R>
  requires std::same_as<std::ranges::range_value_t<R>,
                        std::pair<FastqRecord<>, FastqRecord<>>>
auto
operator|(R&& read_pairs, const BurrowWheelerAligner alinger) {
  // hs37d5 needs to padding some 'A' to the end to avoid SW substr out of
  // bounds since the end of hs37d5 are not 'N's.
  // alinger.ref.seq += istring(alinger.PAIR_DIST, 0);
  auto alignments = std::vector<SamRecord<>>{};
  for (const auto& read_pair : read_pairs) {
    auto [sam1, sam2] = alinger.generate_sam(read_pair);
    alignments.push_back(std::move(sam1));
    alignments.push_back(std::move(sam2));
  }
  return alignments;
}

template<std::ranges::range R>
  requires std::same_as<std::ranges::range_value_t<R>, SamRecord<>>
auto
operator|(R&& alignments, const HaplotypeCaller caller) {
  assert(
    std::ranges::all_of(caller.ref.seq, [](auto c) { return ::isupper(c); }));
  return caller.run(alignments);
}

};  // namespace biovoltron

namespace biovoltron::pipe {

template<int SA_INTV, typename size_type, SASorter Sorter>
using build = FMIndex<SA_INTV, size_type, Sorter>;
using align = BurrowWheelerAligner;
using call = HaplotypeCaller;

};  // namespace biovoltron::pipe
