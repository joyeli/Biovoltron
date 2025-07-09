#pragma once

#include <biovoltron/algo/align/tailor/alignment.hpp>
#include <biovoltron/algo/align/tailor/index.hpp>
#include <biovoltron/file_io/fasta.hpp>
#include <biovoltron/utility/istring.hpp>

namespace biovoltron {

/**
 * todo: some brief descriptions.
 *
 * todo: some detail descriptions.
 *
 * Example:
 * auto sample1 = std::vector<FastqRecord<>>{};
 * {
 *   auto record = FastqRecord{};
 *   auto ifs = std::ifstream{"sample1.fastq"};
 *   while (ifs >> record) sample1.emplace_back(record);
 * }
 *
 * // Load prebuilt index.
 * auto index = Index{};
 * {
 *   auto ifs = std::ifstream{"ref.idx"};
 *   index.load(ifs);
 * }
 *
 * // Load prebuilt index.
 * auto rc_index = Index{};
 * {
 *   auto ifs = std::ifstream{"rc_ref.idx"};
 *   rc_index.load(ifs);
 * }
 *
 * auto tailor = Tailor{index, rc_index};
 * auto alns = std::vector<Alignment>{};
 *
 * for (const auto& read : sample1)
 *   alns.emplace_back(tailor.search(read));
 *
 * auto ofs = std::ofstream{"result.sam"};
 * for (const auto& aln : alns)
 *   for (const auto& record : aln_to_sam_list(aln))
 *     ofs << record << '\n';
 */
template<
  int SA_INTV = 1,
  typename size_type = std::uint32_t,
  SASorter Sorter = PsaisSorter<size_type>
>
struct Tailor {
  using range_type = std::pair<size_type, size_type>;
  using index_type = Index<SA_INTV, size_type, Sorter>;
  constexpr static size_type npos = -1;

  index_type fmi;
  index_type rc_fmi;
  bool allow_seed_mismatch = false;
  bool strict_mode = false;
  bool enable_c2t = false;
  size_type seed_len = 18;
  size_type max_multi = 10;

 private:
  struct Raw {
    bool forward;
    size_type hit_pos;
    range_type rng;
    std::vector<Mismatch> mismatches;
    std::set<std::uint32_t> tc_set;

    auto operator==(const Raw& other) const {
      return hit_pos == other.hit_pos &&
             forward == other.forward &&
             tc_set == other.tc_set &&
             std::ranges::equal(
              mismatches, 
              other.mismatches, 
              {}, 
              &Mismatch::pos, 
              &Mismatch::pos);
    }
  };

  auto exact_match(istring_view read, bool forward) const {
    auto rngs = std::vector<range_type>{};
    const auto& index = forward ? fmi : rc_fmi;

    auto beg = size_type{};
    auto end = index.get_bwt_size();

    while (!read.empty()) {
      if (end <= beg) break;
      beg = index.lf_mapping(read.back(), beg);
      end = index.lf_mapping(read.back(), end);
      rngs.emplace_back(beg, end);
      read.remove_suffix(1);
    }

    auto hit_pos = (end <= beg) ? read.size() + 1 : read.size();
    if (end <= beg) rngs.pop_back();

    return std::tuple(static_cast<size_type>(hit_pos), rngs);
  }

  auto further_scan(
    istring_view read,
    size_type mismatch_pos,
    range_type rng,
    bool forward,
    std::vector<Mismatch> mismatches = {}
  ) const {
      const auto& index = forward ? fmi : rc_fmi;
      auto candidates = std::vector<Raw>{};
      for (auto chr = 0; chr < 4; chr++) {
        if (chr == read[mismatch_pos]) continue;

        auto beg = index.lf_mapping(chr, rng.first);
        auto end = index.lf_mapping(chr, rng.second);
        if (end <= beg) continue;

        auto prev_rng = range_type{};
        auto read_cpy = read;
        read_cpy.remove_suffix(read_cpy.size() - mismatch_pos);
        while (!read_cpy.empty()) {
          if (end <= beg) break;
          else {
            prev_rng.first = beg; 
            prev_rng.second = end;
          }
          beg = index.lf_mapping(read_cpy.back(), beg);
          end = index.lf_mapping(read_cpy.back(), end);
          read_cpy.remove_suffix(1);
        }

        auto hit_pos = (end <= beg) ? read_cpy.size() + 1 : read_cpy.size();
        if (end <= beg) std::tie(beg, end) = prev_rng;

        auto raw = Raw{forward, hit_pos, {beg, end}, {}};
        std::ranges::copy(mismatches, std::back_inserter(raw.mismatches));
        raw.mismatches.emplace_back(mismatch_pos, Codec::to_char(chr));
        candidates.emplace_back(std::move(raw));
      }

      return candidates;
  }

  auto backtrack(
    istring_view read,
    size_type mismatch_pos,
    std::vector<range_type>& match_history,
    bool forward,
    size_type local_seed_len
  ) const {
    auto candidates = std::vector<Raw>{};
    for (; mismatch_pos < read.size(); mismatch_pos++) {
      std::ranges::move(
        further_scan(read, mismatch_pos, match_history.back(), forward),
        std::back_inserter(candidates)
      );
      match_history.pop_back();
    }

    if (candidates.empty())
      return candidates;

    auto min_hit_pos = std::ranges::min_element(
      candidates, {}, &Raw::hit_pos)->hit_pos;
    // If hit position still in seed region, report unmappable.
    if (min_hit_pos >= read.size() - local_seed_len)
      return std::vector<Raw>{};

    if (std::ranges::count(
      candidates, min_hit_pos, &Raw::hit_pos) > max_multi)
      return std::vector<Raw>{};
      
    auto final_candidates = decltype(candidates){};
    std::ranges::move(
      candidates | std::views::filter([min_hit_pos](const Raw& raw) { 
        return raw.hit_pos == min_hit_pos; }), 
      std::back_inserter(final_candidates)
    );
    return final_candidates;
  }

  auto search_(istring_view read, bool forward, size_type local_seed_len) const {
    auto [hit_pos, match_history] = exact_match(read, forward);

    auto candidates = std::vector<Raw>{};
    if (hit_pos < read.size() - local_seed_len)
      candidates.emplace_back(forward, hit_pos, match_history.back());
    else if (allow_seed_mismatch)
      std::ranges::move(
        backtrack(read, hit_pos-1, match_history, forward, local_seed_len),
        std::back_inserter(candidates)
      );

    if (candidates.empty() || candidates.front().hit_pos <= 1)
      return candidates;

    auto final_candidates = std::vector<Raw>{};
    for (auto& candidate : candidates) {
      auto results = further_scan(
        read, 
        candidate.hit_pos-1, // Mismatch pos.
        candidate.rng, 
        forward, 
        candidate.mismatches);

      auto hit_zero = [](size_type hit_pos) { return hit_pos == 0; };
      if (std::ranges::none_of(results, hit_zero, &Raw::hit_pos))
        // Found tail.
        final_candidates.emplace_back(
          forward, candidate.hit_pos, candidate.rng, candidate.mismatches);
      else
        std::ranges::move(results, std::back_inserter(final_candidates));
    }

    return final_candidates;
  }

  auto pick_best(std::vector<Raw>& candidates) const {

    auto got_best = [](const auto& candidates) {
      return std::equal(candidates.cbegin() + 1, candidates.cend(), candidates.cbegin());
    };

    if (candidates.empty() || got_best(candidates))
      return candidates;

    // Pick shortest tail.
    {
      const auto min_hit_pos = std::ranges::min(candidates, {}, &Raw::hit_pos).hit_pos;
      const auto ret = std::ranges::remove_if(candidates, 
        [min_hit_pos](const auto& raw) { return raw.hit_pos != min_hit_pos; });
      candidates.erase(ret.begin(), ret.end());
    }

    // Pick fewest mismatches.
    if (!got_best(candidates)) {
      const auto min_mismatch_num = std::ranges::min(
        candidates, 
        [](const auto& lhs, const auto& rhs) { return lhs.size() < rhs.size(); },
        &Raw::mismatches).mismatches.size();
      const auto ret = std::ranges::remove_if(candidates, 
        [min_mismatch_num](const auto& raw) { return raw.mismatches.size() != min_mismatch_num; });
      candidates.erase(ret.begin(), ret.end());
    }
                
    // Pick mismatch toward 3' (Note: this can be discussed).
    if (!got_best(candidates) && !candidates.front().mismatches.empty()) {
      const auto mm = std::ranges::min(
        candidates, 
        [](const auto& lhs, const auto& rhs) { 
          return std::ranges::max(lhs, {}, &Mismatch::pos).pos < std::ranges::max(rhs, {}, &Mismatch::pos).pos; },
        &Raw::mismatches).mismatches;
      const auto smallest_max_mismatch_pos = std::ranges::max(mm, {}, &Mismatch::pos).pos;
      const auto ret = std::ranges::remove_if(candidates, 
        [smallest_max_mismatch_pos](const auto& raw) { 
          return std::ranges::max(raw.mismatches, {}, &Mismatch::pos).pos != smallest_max_mismatch_pos; });
      candidates.erase(ret.begin(), ret.end());
    }

    // Note: this filter doesn't make sense
    // // Prefer forward strand(Note: this can be discussed).
    // if (!got_best(candidates)) {
    //   // Forward stand is "false" because we did a reverse complement.
    //   // So we remove "true".
    //   const auto ret = std::ranges::remove(candidates, true, &Raw::forward);
    //   candidates.erase(ret.begin(), ret.end());
    // }

    return candidates;
  }

  template <bool Encoded>
  auto raws2alignment(const FastqRecord<Encoded>& record, std::vector<Raw>&& raws) const {
    auto aln = Alignment{};
    aln.name = record.name;
    aln.seq = record.seq;
    aln.qual = record.qual;

    // Report unmappable with empty hits.
    if (raws.empty()) 
      return aln;

    const auto forward = raws.front().forward;
    const auto hit_pos = raws.front().hit_pos;
    const auto& index = forward ? fmi : rc_fmi;
    const auto& read_len = record.seq.size();
    auto get_reverse = [](size_type pos, size_type total_len) {
      return total_len - pos - 1;
    };

    // Reverse the result back.
    aln.forward = !forward;
    aln.tail_pos = (hit_pos == 0 ) ? -1 : get_reverse(hit_pos-1, read_len);
    const auto interval_len = (aln.tail_pos == -1) ? read_len : aln.tail_pos;

    for (auto& raw : raws) {
      const auto [beg, end] = raw.rng;
      for (auto& iv : index.get_intervals(beg, end, interval_len)) {
        // Reverse the result back.
        if (aln.forward) {
          iv.begin--; iv.end--; std::swap(iv.begin, iv.end);
          iv.begin = get_reverse(iv.begin, index.get_chr_size(iv.chrom));
          iv.end = get_reverse(iv.end, index.get_chr_size(iv.chrom));
        }
        iv.strand = aln.forward ? '+' : '-';

        // Reverse mismatch positions and reverse complement correct base.
        for (auto& mm : raw.mismatches) {
          mm.correct_base = Codec::comp(mm.correct_base);
          mm.pos = get_reverse(mm.pos, read_len);
        }

        aln.hits.emplace_back(raw.mismatches, raw.tc_set, iv);
      }
    }
    return aln;
  }

  template <bool Encoded>
  auto c2t (
    const FastqRecord<Encoded>& record, 
    std::vector<Raw>& candidates,
    size_type local_seed_len
  ) const {
    // todo: Make combination cleaner and more expressive.
    auto read = record.seq;
    auto pos = std::vector<size_type>{};
    for (auto i = 0; i < read.size(); i++)
      if (read[i] == 'C')
        pos.push_back(i);
    // At most do 3-combination
    for (auto i = 1; i <= 3 && i <= pos.size(); i++) {
    
    // exhausted all combination
    // for (auto i = 1; i <= pos.size(); i++) {
      auto list = std::string(i, 1);
      list.resize(pos.size(), 0);
  
      do {
        auto indexes = std::vector<size_type>{};
        for (auto j = 0; j < list.size(); j++)
          if (list[j] == 1)
            indexes.push_back(pos[j]);
  
        for (auto idx : indexes)
          read[idx] = 'T';
        auto rc_read = istring{};
        if constexpr (Encoded) rc_read = Codec::rev_comp(read);
        else rc_read = Codec::rev_comp(Codec::to_istring(read));

        std::ranges::move(search_(rc_read, true, local_seed_len), std::back_inserter(candidates));
        std::ranges::move(search_(rc_read, false, local_seed_len), std::back_inserter(candidates));
  
        if (!candidates.empty()) {
          for (auto& candidate : candidates) {
            candidate.tc_set = std::set<std::uint32_t>(indexes.begin(), indexes.end());
          }
          return;
        }
  
        for (auto idx : indexes)
          read[idx] = 'C';
      } while (std::ranges::prev_permutation(list).found);
    }
  }

 public:
  Tailor(): fmi(), rc_fmi() {}
  
  Tailor(const index_type& index, const index_type& rc_index) 
    : fmi(std::move(index)), rc_fmi(std::move(rc_index)) {
    assert(fmi.get_bwt_size() != 0);
    assert(rc_fmi.get_bwt_size() != 0);
  }

  template <bool Encoded>
  auto search(const FastqRecord<Encoded>& record) const {
    size_type local_seed_len = seed_len;

    auto make_fallback_return = [](bool is_N){
      return std::make_pair(
        Alignment{.seq = (is_N ? "N": "")},
        Alignment{.seq = (is_N ? "N": "")}
      );
    };

    if (record.seq.size() < local_seed_len)
      return make_fallback_return(false);

    auto rc_read = istring{};
    if constexpr (Encoded) rc_read = Codec::rev_comp(record.seq);
    else rc_read = Codec::rev_comp(Codec::to_istring(record.seq));

    // make aln of reads that contain 'N' a lil bit special
    // so the aligner can calculate the filtered read count
    if (std::ranges::find(rc_read, Codec::to_int('N')) != rc_read.end())
      return make_fallback_return(true);

    auto candidates = std::vector<Raw>{};
    std::ranges::move(search_(rc_read, true, local_seed_len), std::back_inserter(candidates));
    std::ranges::move(search_(rc_read, false, local_seed_len), std::back_inserter(candidates));

    if (candidates.empty() && enable_c2t) 
      c2t(record, candidates, local_seed_len);

    if (!candidates.empty() && strict_mode) {
      // todo:Apply strict mode.
      
      auto min_hit_pos = std::ranges::min(candidates, {}, &Raw::hit_pos).hit_pos;
      if(min_hit_pos > 5) { // still no tail smaller than 5
        local_seed_len = rc_read.size() - min_hit_pos;

        candidates.clear();

        std::ranges::move(search_(rc_read, true, local_seed_len), std::back_inserter(candidates));
        std::ranges::move(search_(rc_read, false, local_seed_len), std::back_inserter(candidates));

        if (candidates.empty() && enable_c2t)
          c2t(record, candidates, local_seed_len);
      }
    }
    
    auto&& best = pick_best(candidates);
    auto filter_strand = [](const std::vector<Raw>& candidates, bool forward){
      std::vector<Raw> filtered_candidates;
      std::ranges::copy_if(
        candidates,
        std::back_inserter(filtered_candidates),
        [&forward](bool f){ return forward == f; },
        &Raw::forward
      );
      return filtered_candidates;
    };
    
    // Note that forward == false means read is aln to forward strand
    auto aln_forward = raws2alignment(record, filter_strand(best, false));
    auto aln_reverse = raws2alignment(record, filter_strand(best, true));
    
    if ((aln_forward.hits.size() + aln_reverse.hits.size()) > max_multi)
      return make_fallback_return(false);

    return std::make_pair(aln_forward, aln_reverse);
  }
};

}  // namespace biovoltron
