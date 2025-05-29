#pragma once

#include <biovoltron/algo/align/exact_match/fm_index.hpp>
#include <biovoltron/algo/align/inexact_match/smithwaterman_sse.hpp>
#include <biovoltron/algo/align/mapq/mapq.hpp>
#include <biovoltron/algo/sort/stable_sorter.hpp>
#include <biovoltron/file_io/fasta.hpp>
#include <biovoltron/file_io/fastq.hpp>
#include <biovoltron/file_io/sam.hpp>
#include <map>
#include <spdlog/spdlog.h>
#include <sstream>

// Note: Be aware, this aligner is specifically for
// hs37d5 dataset. It is not recommended to use it for
// other dataset.

namespace biovoltron {

/**
 * @ingroup applications
 * A highly optimized Aligner (implement using modern C++20) 
 * specifically for paired-end short read alignment with read 
 * length ~150 bp.
 *
 * **This aligner is highly optimized on the following sequencing 
 * characteristic (other datasets are not recommended)**:
 * - read length: **148~150bp**
 * 
 * If you know the mean and variance of the insert size of the 
 * sequencing data, we highly recommend you pass it into the aligner.
 */
struct BurrowWheelerAligner {
  struct Parameters {
    const int INSERT_MEAN = 550;
    const int INSERT_VAR = 150;
    const int PAIR_DIST = INSERT_MEAN + 4 * INSERT_VAR + 50;

    const int MAX_HIT_CNT = 512;
    const int MAX_EM_CNT = 128;
    const int MAX_SW_CNT = 32;
    const int MAX_RESCUE_CNT = 128;
    const int MAX_SEED_CNT = 4;

    const int SEED_LEN = 19;
    const int SEED_OVERLAP = 4;
    const int EXTEND = 100;
    const int SW_THRESHOLD = 30;
    const int KMER_SIZE = 8;
    const int MIN_FIND_CNT = 4;
    const int MAX_FIND_CNT_DIFF = 4;
    const int MAX_SW_DIFF = 30;
    const int PEN_UNPAIRED = 19;
  };
  const FastaRecord<true> ref;
  const FMIndex<1, uint32_t, StableSorter<uint32_t>> index;
  const Parameters args;

  auto
  print() {
    SPDLOG_DEBUG("================== argument ==================");
    SPDLOG_DEBUG("MAX_HIT_CNT: {}", args.MAX_HIT_CNT);
    SPDLOG_DEBUG("MAX_EM_CNT: {}", args.MAX_EM_CNT);
    SPDLOG_DEBUG("MAX_SW_CNT: {}", args.MAX_SW_CNT);
    SPDLOG_DEBUG("MAX_RESCUE_CNT: {}", args.MAX_RESCUE_CNT);
    SPDLOG_DEBUG("MAX_SEED_CNT: {}", args.MAX_SEED_CNT);

    SPDLOG_DEBUG("SEED_LEN: {}", args.SEED_LEN);
    SPDLOG_DEBUG("SEED_OVERLAP: {}", args.SEED_OVERLAP);
    SPDLOG_DEBUG("EXTEND: {}", args.EXTEND);
    SPDLOG_DEBUG("SW_THRESHOLD: {}", args.SW_THRESHOLD);
    SPDLOG_DEBUG("KMER_SIZE: {}", args.KMER_SIZE);
    SPDLOG_DEBUG("MIN_FIND_CNT: {}", args.MIN_FIND_CNT);
    SPDLOG_DEBUG("MAX_FIND_CNT_DIFF: {}", args.MAX_FIND_CNT_DIFF);
    SPDLOG_DEBUG("MAX_SW_DIFF: {}", args.MAX_SW_DIFF);
    SPDLOG_DEBUG("PEN_UNPAIRED: {}", args.PEN_UNPAIRED);
    SPDLOG_DEBUG("INSERT_MEAN: {}", args.INSERT_MEAN);
    SPDLOG_DEBUG("INSERT_VAR: {}", args.INSERT_VAR);
    SPDLOG_DEBUG("PAIR_DIST: {}", args.PAIR_DIST);
  }

  struct Aln {
    std::uint32_t pos{};
    std::uint8_t score{};
    std::uint8_t score2{};
    bool forward{1};
    std::uint8_t read_end{};
    std::uint32_t ref_end{};
    std::uint8_t find_cnt{};
    std::uint8_t align_len{};
    std::uint8_t mapq{};
    std::uint8_t sub_score{};
    bool rescued{};
    std::string cigar;
    istring rev_comp;

    friend auto&
    operator<<(std::ostream& os, Aln aln) {
      return os << "(pos: " << aln.pos << (aln.forward ? "(->)" : "(<-)")
                << ", score: " << +aln.score << ", score2: "
                << +aln.score2
                // << ", align_len: " << +aln.align_len << ", find_cnt: " <<
                // +aln.find_cnt
                << ", rescued:" << aln.rescued << ", mapq:" << +aln.mapq << ")";
    }

    auto
    operator==(const Aln& other) const {
      return pos == other.pos && forward == other.forward;
    }

    auto
    operator<=>(const Aln& other) const {
      if (auto cmp = pos <=> other.pos; cmp != 0)
        return cmp;
      if (auto cmp = forward <=> other.forward; cmp != 0)
        return cmp;
      if (auto cmp = other.score <=> score; cmp != 0)
        return cmp;
      return other.cigar.size() <=> cigar.size();
    }
  };

  struct Anchor {
    std::uint32_t ref_pos{};
    std::uint8_t seed_pos{};
    std::uint8_t seed_size{};
    bool forward{};
    bool repeat{};

    auto
    operator<=>(const Anchor& other) const noexcept = default;
    friend auto&
    operator<<(std::ostream& os, const Anchor& an) {
      return os << "(" << +an.seed_pos << ": " << an.ref_pos << ", "
                << +an.seed_size << ")";
    }
  };

  struct AlnPair {
    Aln aln1;
    Aln aln2;
    auto
    dist() const noexcept {
      return DIFF(aln1.pos, aln2.pos);
    }
    auto
    score() const noexcept {
      return aln1.score + aln2.score;
    }
  };

  struct SeedSpan {
    istring_view seed;
    std::span<const std::uint32_t> span;

    auto
    operator<=>(const SeedSpan& other) const noexcept {
      if (auto cmp = span.size() <=> other.span.size(); cmp != nullptr)
        return cmp;
      return other.seed.size() <=> seed.size();
    }
    bool
    operator==(const SeedSpan& other) const noexcept {
      return seed.size() == other.seed.size()
             && span.size() == other.span.size();
    }
  };

 private:
  constexpr static auto DIFF
    = [](auto a, auto b) { return (a > b) ? (a - b) : (b - a); };

  auto
  filter_alns(std::vector<Aln>& alns) const {
    const auto best_score = alns.front().score;
    const auto last = std::ranges::find_if(
      alns,
      [best_score, MAX_SW_DIFF = this->args.MAX_SW_DIFF](const auto score) {
        return score < best_score - MAX_SW_DIFF;
      },
      &Aln::score);
    const auto final_size = std::ranges::distance(alns.begin(), last);
    alns.resize(final_size);

    for (auto i = 0; const auto& aln : alns | std::views::drop(final_size)) {
      SPDLOG_DEBUG("\n************ filtered results ************");
      SPDLOG_DEBUG("[{}] {}", i++, aln);
    }
  }

  auto
  finalize_alns(std::vector<Aln>& alns) const {
    if (alns.size() <= 1)
      return;
    std::ranges::sort(alns);
    const auto [begin, end] = std::ranges::unique(alns);
    alns.erase(begin, end);
    std::ranges::sort(alns, std::ranges::greater{}, &Aln::score);
    filter_alns(alns);
  }

  static auto
  print_alns(const std::vector<Aln>& alns) {
    for (auto i = 0; const auto& aln : alns | std::views::take(32))
      SPDLOG_DEBUG("[{}] {}", i++, aln);
  }

  auto
  split_read(istring_view read) const {
    constexpr auto npos = istring_view::npos;
    auto frags = std::vector<istring_view>{};
    for (auto start = npos, i = 0ul; i <= read.size(); i++) {
      if (i == read.size() || read[i] == 4) {
        if (start != npos && i - start >= args.SEED_LEN)
          frags.push_back(read.substr(start, i - start));
        start = npos;
      } else if (start == npos)
        start = i;
    }
    return frags;
  }

  static auto
  need_revert(istring_view read, istring_view ref) {
    const auto max_mis_cnt = (read.size() + 4) / 5;
    auto mis_cnt = 0;
    for (auto i = 0; mis_cnt <= max_mis_cnt && i < read.size(); i++)
      if (read[i] != ref[i] && read[i] != 4 && ref[i] != 4)
        mis_cnt++;
    return mis_cnt <= max_mis_cnt;
  }

  auto
  set_cigar(Aln& aln, istring_view read, auto& profile) const {
    SPDLOG_DEBUG("============== compute cigar ==============");

    if (!aln.cigar.empty())
      return;
    if (aln.score == read.size()) {
      aln.cigar = std::to_string(read.size()) + 'M';
      SPDLOG_DEBUG("full score: {}", aln.cigar);
      aln.align_len = read.size();
      return;
    }

    const auto sw_pos = aln.ref_end - read.size() - args.EXTEND;
    const auto subref
      = istring_view{ref.seq}.substr(sw_pos, read.size() + args.EXTEND + 1);

    if (profile.read == nullptr)
      profile = SseSmithWaterman::get_profile(read);
    auto sw
      = SseSmithWaterman::align(profile, subref, true, true, args.SW_THRESHOLD);

    const auto ref_beg = sw_pos + sw.ref_beg;
    const auto ref_end = sw_pos + sw.ref_end;
    aln.pos = ref_beg;
    aln.score = sw.score;
    aln.cigar = std::move(sw.cigar);
    aln.align_len
      = std::max(sw.ref_end - sw.ref_beg + 1, sw.read_end - sw.read_beg + 1);

    SPDLOG_DEBUG("pos: {}", aln.pos);
    SPDLOG_DEBUG("raw cigar: {}", aln.cigar);
  }

  auto
  display_chains(const auto& chains) const {
    SPDLOG_DEBUG("************ seed chain ************");
	std::stringstream ss;
    for (const auto& anchors : chains | std::views::take(8)) {
      const auto& front = anchors.front();
      ss << "(" << anchors.size() << ", "
                << (front.forward ? "->): " : "<-): ");
      ss << front;
      for (const auto anchor : anchors | std::views::drop(1))
        ss << "->" << anchor;
      SPDLOG_DEBUG("{}", ss.str());
      ss.str("");
    }
  }

  static auto
  get_score(istring_view read, istring_view ref, bool forward)
    -> std::pair<std::uint16_t, Cigar> {
    const auto full_score = read.size();

    const auto read_beg = read.substr(0, 5);
    const auto ref_beg = ref.substr(0, 5);
    const auto read_end = read.substr(read.size() - 5);
    const auto ref_end = ref.substr(ref.size() - 5);

    read.remove_prefix(5);
    ref.remove_prefix(5);
    read.remove_suffix(5);
    ref.remove_suffix(5);

    constexpr auto only_one_mismatch = [](istring_view read, istring_view ref) {
      const auto [in1, in2] = std::ranges::mismatch(read, ref);
      return read.substr(in1 - read.begin() + 1)
             == ref.substr(in2 - ref.begin() + 1);
    };

    if (read != ref) {
      if (read_beg != ref_beg || read_end != ref_end)
        return {};
      if (!only_one_mismatch(read, ref))
        return {};

      return {full_score - 5, std::to_string(full_score) + 'M'};
    }

    if (read_beg == ref_beg && read_end == ref_end) {
      return {full_score, std::to_string(full_score) + 'M'};
    }

    constexpr auto get_beg_score =
      [](auto read_beg, auto ref_beg) -> std::pair<std::uint16_t, std::string> {
      constexpr static auto cigars
        = std::array{"5S", "4S1M", "3S2M", "2S3M", "1S4M"};
      for (auto clip_idx = 4; clip_idx >= 0; clip_idx--) {
        if (read_beg[clip_idx] != ref_beg[clip_idx]) {
          auto cigar = cigars[4 - clip_idx];

          return {4 - clip_idx, cigar};
        }
      }
      return {5, "5M"};
    };

    constexpr auto get_end_score =
      [](auto read_end, auto ref_end) -> std::pair<std::uint16_t, std::string> {
      constexpr static auto cigars
        = std::array{"5S", "1M4S", "2M3S", "3M2S", "4M1S"};
      for (auto clip_idx = 0; clip_idx <= 4; clip_idx++) {
        if (read_end[clip_idx] != ref_end[clip_idx]) {
          auto cigar = cigars[clip_idx];
          return {clip_idx, cigar};
        }
      }
      return {5, "5M"};
    };

    auto score = read.size();
    auto mid_cigar = std::to_string(score) + 'M';
    const auto [end_score, end_cigar] = get_end_score(read_end, ref_end);
    const auto [beg_score, beg_cigar] = get_beg_score(read_beg, ref_beg);
    score += (beg_score + end_score);
    if (score < full_score - 5)
      return {};
    auto cigar = static_cast<Cigar>(beg_cigar + mid_cigar + end_cigar);
    cigar.compact();
    return {score, std::move(cigar)};
  }

  static auto
  compute_se_mapq(const std::vector<Aln>& alns, float frac_rep) {
    const auto [opt_score, sub_score, sub_cnt]
      = get_opt_subopt_count(alns | std::views::transform(&Aln::score));
    const auto& front = alns.front();
    const auto mapq = mem_approx_mapq_se(
      {opt_score, front.score2, sub_score, front.align_len, sub_cnt, frac_rep});
    return std::pair{mapq, sub_score};
  }

  auto
  get_best_one(const std::vector<Aln>& alns, istring_view read,
               istring_view rread, auto& profile, auto& rprofile,
               float frac_rep) const {
    auto aln = alns.front();
    if (aln.forward)
      set_cigar(aln, read, profile);
    else
      set_cigar(aln, rread, rprofile);
    const auto [mapq, sub_score] = compute_se_mapq(alns, frac_rep);
    aln.mapq = mapq;
    aln.sub_score = sub_score;
    return aln;
  }

  // O(NlogN) pairing
  auto
  pairing2(std::vector<Aln> alns1, std::vector<Aln> alns2) const {
    std::ranges::sort(alns1);
    std::ranges::sort(alns2);
    auto aln_pairs = std::vector<AlnPair>{};
    auto begin = alns2.begin();
    auto end = begin;
    for (const auto& aln1 : alns1) {
      while (begin != alns2.end() && begin->pos < aln1.pos - args.PAIR_DIST)
        ++begin;
      while (end != alns2.end() && end->pos < aln1.pos + args.PAIR_DIST) ++end;
      for (auto it = begin; it != end; ++it) {
        const auto& aln2 = *it;
        if (aln1.forward != aln2.forward)
          aln_pairs.emplace_back(aln1, aln2);
      }
    }
    std::ranges::sort(aln_pairs, std::ranges::greater{}, &AlnPair::score);
    return aln_pairs;
  }

  static auto
  print_paires(const std::vector<AlnPair>& aln_pairs) {
    for (auto i = 0; const auto& aln_pair : aln_pairs) {
      SPDLOG_DEBUG("[{}]{{ ({}) <-- {} --> ({}) }} (score: {})", i++, aln_pair.aln1, aln_pair.dist()
        , aln_pair.aln2, aln_pair.score());
    }
  }

  auto
  get_equal_span(std::span<const std::uint32_t> span, std::uint32_t offset,
                 istring_view equal_seed) const {
    const auto [begin, end] = std::ranges::equal_range(
      span, equal_seed, {},
      [this, offset, equal_size = equal_seed.size()](const auto pos) {
        return istring_view{ref.seq}.substr(pos + offset, equal_size);
      });
    return span.subspan(begin - span.begin(), end - begin);
  }

  auto
  get_spans(istring_view read) const {
    auto seed_spans = std::vector<SeedSpan>{};
    auto repeat_size = 0;
    auto origin_read = read;
    while (read.size() >= args.SEED_LEN) {
      const auto seed = read.substr(read.size() - args.SEED_LEN);
      const auto [begin, end, offset] = index.get_range(seed, 0);
      const auto span = index.get_offsets(begin, end);
      if (span.size() <= args.MAX_HIT_CNT) {
        const auto cur_seed = seed.substr(offset);
        SPDLOG_DEBUG("seed: {}({}) -> ({})", cur_seed, cur_seed.size(), span.size());
        if (!span.empty())
          seed_spans.emplace_back(seed.substr(offset), span);
        read.remove_suffix(args.SEED_LEN - offset - args.SEED_OVERLAP);
      } else {
        assert(offset == 0);
        const auto seed2 = read.substr(0, read.size() - args.SEED_LEN);
        const auto [begin2, end2, offset2]
          = index.get_range(seed2, begin, end, args.MAX_HIT_CNT);
        const auto span2 = index.get_offsets(begin2, end2);
        if (span2.size() <= args.MAX_HIT_CNT) {
          const auto cur_seed = read.substr(offset2);
          SPDLOG_DEBUG("seed: {}({}) -> ({})", cur_seed, cur_seed.size(), span2.size());
          if (!span2.empty())
            seed_spans.emplace_back(read.substr(offset2), span2);
          read = read.substr(0, offset2 + args.SEED_OVERLAP);
        } else {
          assert(offset2 == 0);
          SPDLOG_DEBUG("seed: {}({}) -> ({})". read, read.size(), span2.size());
          SPDLOG_DEBUG("\n************ backtrace extend ************");

          auto remain_seed = origin_read.substr(read.size());
          auto equal_offset = read.size();
          auto equal_span = span2;
          while (!remain_seed.empty()) {
            const auto equal_size
              = std::min(static_cast<int>(remain_seed.size()),
                         static_cast<int>(std::log2(equal_span.size()) / 2));
            const auto equal_seed = remain_seed.substr(0, equal_size);
            equal_span = get_equal_span(equal_span, equal_offset, equal_seed);

            SPDLOG_DEBUG("equal_seed: {}({}) -> ({})", equal_seed, equal_seed.size(), equal_span.size());
            remain_seed.remove_prefix(equal_size);
            equal_offset += equal_size;
            if (equal_span.size() <= args.MAX_HIT_CNT) {
              const auto extend_seed = origin_read.substr(0, equal_offset);
              seed_spans.emplace_back(extend_seed, equal_span);

              SPDLOG_DEBUG("seed: {}({}) -> ({})", extend_seed, extend_seed.size(), equal_span.size());

              break;
            }
          }

          repeat_size = read.size();
          read = {};
          break;
        }
      }
    }
    // recover
    if (read.size() >= args.SEED_LEN - args.SEED_OVERLAP) {
      const auto [begin, end, offset] = index.get_range(read, 0);
      const auto span = index.get_offsets(begin, end);
      if (span.size() <= args.MAX_HIT_CNT) {
        SPDLOG_DEBUG("seed: {} -> ({})", read.substr(offset), span.size());
        if (!span.empty())
          seed_spans.emplace_back(read.substr(offset), span);
      }
    }
    return std::pair{std::move(seed_spans), repeat_size};
  }

  auto
  seeding_impl(istring_view read, bool forward) const {
    const auto frags = split_read(read);
    SPDLOG_DEBUG("\n************ read seeds ************");

    auto seed_spans = std::vector<SeedSpan>{};
    auto repeats = 0;
    for (const auto frag : frags) {
      const auto [spans, repeat_size] = get_spans(frag);
      repeats += repeat_size;
      for (const auto span : spans) seed_spans.push_back(span);
    }
    std::ranges::sort(seed_spans);

    auto anchors = std::vector<Anchor>{};
    for (const auto [seed, span] :
         seed_spans | std::views::take(args.MAX_SEED_CNT)) {
      const auto seed_pos = seed.data() - read.data();
      for (const auto ref_pos : span)
        anchors.emplace_back(ref_pos, seed_pos, seed.size(), forward);
    }

    auto chain_map = std::map<std::uint32_t, std::vector<Anchor>>{};
    for (const auto anchor : anchors) {
      const auto read_pos = anchor.ref_pos - anchor.seed_pos;
      const auto lower = chain_map.lower_bound(read_pos - args.SEED_LEN);
      const auto upper = chain_map.upper_bound(read_pos + args.SEED_LEN);
      if (lower != upper) {
        for (auto it = lower; it != upper; ++it) it->second.push_back(anchor);
        continue;
      }
      chain_map[read_pos].push_back(anchor);
    }

    auto chains = std::vector<std::vector<Anchor>>{};
    for (auto& chain : chain_map | std::views::values) {
      std::ranges::sort(chain);
      chains.push_back(std::move(chain));
    }
    return std::pair{std::move(chains), repeats};
  }

  auto
  seeding(istring_view read, istring_view rread) const {
    auto [chains, repeats] = seeding_impl(read, true);

    SPDLOG_DEBUG("\n************ reverse ************");
    std::stringstream ss;
    ss << "rread: ";
    for (auto i = 0; i < rread.size(); i++)
      ss << (rread[i] != 4 ? Codec::to_char(rread[i]) : '|');
    SPDLOG_DEBUG("{}", ss.str());

    const auto [rchains, rrepeats] = seeding_impl(rread, false);
    std::ranges::copy(rchains, std::back_inserter(chains));
    std::ranges::sort(chains, std::ranges::greater{},
                      &std::vector<Anchor>::size);
    return std::tuple{std::move(chains), repeats, rrepeats};
  }

  auto
  exact_match(auto& chains, istring_view read, istring_view rread,
              int find_cnt) const {
    auto alns = std::vector<Aln>{};
    auto sw_chains = std::vector<std::vector<Anchor>>{};
    for (auto& chain : chains) {
      const auto [ref_pos, seed_pos, seed_size, forward, repeat]
        = chain.front();

      auto read_pos = ref_pos - seed_pos;
      const auto sw_read = forward ? read : rread;

      SPDLOG_DEBUG("ref_pos: {}", ref_pos);
      SPDLOG_DEBUG("read_pos: {}", read_pos);
      SPDLOG_DEBUG("sw_read: {}", sw_read);
      if (const auto [score, cigar] = get_score(
            sw_read, istring_view{ref.seq}.substr(read_pos, read.size()),
            forward);
          score) {
        if (const auto [size, op] = cigar.front(); op == 'S') {
          read_pos += size;
        }
        alns.emplace_back(read_pos, score, 0, forward).cigar = cigar;
        alns.back().find_cnt = find_cnt;
        alns.back().align_len = cigar.ref_size();
      } else {
        sw_chains.push_back(std::move(chain));
      }
    }
    if (sw_chains.size() > args.MAX_EM_CNT)
      sw_chains.resize(args.MAX_EM_CNT);
    return std::pair{std::move(alns), std::move(sw_chains)};
  }

  auto
  get_kmers(istring_view read) const {
    auto kmers = std::vector<std::uint32_t>{};
    // overlap one base
    for (auto i = 0; i < read.size() / (args.KMER_SIZE - 1); i++) {
      auto needle = read.substr(i * (args.KMER_SIZE - 1), args.KMER_SIZE);
      kmers.push_back(Codec::hash(needle));
    }
    return kmers;
  }

  auto
  find_kmers(const std::vector<std::uint32_t>& kmers, istring_view ref,
             std::vector<bool>& table) const {
    table.clear();
    table.resize(1 << args.KMER_SIZE * 2);

    for (auto i = 0; i < ref.size() - args.KMER_SIZE + 1; i++)
      table[Codec::hash(ref.substr(i, args.KMER_SIZE))] = true;

    auto find_cnt = 0;
    for (const auto kmer : kmers) find_cnt += table[kmer];
    return find_cnt;
  }

  auto
  get_sw_alns(const auto& chains, int read_size,
              const std::vector<std::uint32_t>& kmers,
              const std::vector<std::uint32_t>& rkmers,
              std::vector<bool>& table, int min_find_cnt) const {
    auto alns = std::vector<Aln>{};
    for (const auto& chain : chains) {
      const auto [ref_pos, seed_pos, seed_size, forward, repeat]
        = chain.front();
      const auto read_pos = ref_pos - seed_pos;
      const auto front_pad
        = seed_pos <= args.EXTEND / 2 ? seed_pos * 2 : args.EXTEND;
      const auto sw_pos = read_pos - front_pad;
      const auto subref
        = istring_view{ref.seq}.substr(sw_pos, read_size + 2 * args.EXTEND);

      const auto find_cnt = find_kmers(forward ? kmers : rkmers, subref, table);
      if (find_cnt < min_find_cnt)
        continue;
      min_find_cnt = std::max(find_cnt - args.MAX_FIND_CNT_DIFF, min_find_cnt);

      alns.emplace_back(sw_pos, 0, 0, forward, 0, 0, find_cnt);
    }
    std::ranges::sort(alns, std::ranges::greater{}, &Aln::find_cnt);
    return std::pair{std::move(alns), min_find_cnt};
  }

  auto
  get_sw_candidates(bool alns_empty, const auto& chains, int read_size,
                    const std::vector<std::uint32_t>& kmers,
                    const std::vector<std::uint32_t>& rkmers,
                    std::vector<bool>& table) const {
    if (alns_empty) {
      return get_sw_alns(chains, read_size, kmers, rkmers, table,
                         args.MIN_FIND_CNT);
    } else {
      const auto indel_chains
        = chains | std::views::take_while([this](const auto& chain) {
            return chain.size() >= args.MAX_SEED_CNT / 2;
          });
      return get_sw_alns(indel_chains, read_size, kmers, rkmers, table,
                         kmers.size() - args.MAX_FIND_CNT_DIFF);
    }
  }

  auto
  extending(std::vector<Aln>& alns, const std::vector<Aln>& sw_alns,
            istring_view read, istring_view rread) const {
    auto profile = s_profile{}, rprofile = s_profile{};
    for (auto min_score = args.SW_THRESHOLD; const auto& sw_aln : sw_alns) {
      const auto subref = istring_view{ref.seq}.substr(
        sw_aln.pos, read.size() + 2 * args.EXTEND);
      auto sw = SseSmithWaterman::SWResult{};
      if (sw_aln.forward) {
        if (profile.read == nullptr)
          profile = SseSmithWaterman::get_profile(read);
        sw = SseSmithWaterman::align(profile, subref, false, false, min_score);
      } else {
        if (rprofile.read == nullptr)
          rprofile = SseSmithWaterman::get_profile(rread);
        sw = SseSmithWaterman::align(rprofile, subref, false, false, min_score);
      }

      const auto score = static_cast<int>(sw.score);
      if (score < min_score)
        continue;
      const auto ref_end = sw_aln.pos + sw.ref_end;
      alns.emplace_back(ref_end - sw.read_end, score, sw.score2, sw_aln.forward,
                        sw.read_end, ref_end, sw_aln.find_cnt);

      min_score = std::max(min_score, score - args.MAX_SW_DIFF);
    }
    return std::pair{std::move(profile), std::move(rprofile)};
  }

  static auto
  release_memory(auto& v) {
    v.clear();
    v.shrink_to_fit();
  }

  auto
  rescue(const std::vector<Aln>& alns1, const std::vector<Aln>& alns2,
         istring_view read2, istring_view rread2, auto& profile2,
         auto& rprofile2, const std::vector<std::uint32_t>& kmers2,
         const std::vector<std::uint32_t>& rkmers2, std::vector<bool>& table,
         int min_find_cnt) const {
    const auto [opt_score, sub_score, sub_cnt]
      = get_opt_subopt_count(alns1 | std::views::transform(&Aln::score));
    const auto rescue_cnt = std::min(sub_cnt + 1, args.MAX_RESCUE_CNT);

    SPDLOG_DEBUG("============== rescue count: {} ==============", rescue_cnt);

    auto rescues = std::vector<Aln>{};
    for (auto min_score = args.SW_THRESHOLD;
         const auto& aln1 : alns1 | std::views::take(rescue_cnt)) {
      const auto pos1 = aln1.pos;
      if (std::ranges::any_of(
            alns2,
            [pos1, PAIR_DIST = this->args.PAIR_DIST](const auto pos2) {
              return DIFF(pos1, pos2) <= PAIR_DIST;
            },
            &Aln::pos)) {
        SPDLOG_DEBUG("pos: {} already seen.", pos1);
        continue;
      }

      const auto forward1 = aln1.forward;
      const auto sw_pos = forward1 ? pos1 - args.EXTEND : pos1 - args.PAIR_DIST;
      const auto subref = istring_view{ref.seq}.substr(
        sw_pos, args.EXTEND + read2.size() + args.PAIR_DIST);

      const auto find_cnt
        = find_kmers(forward1 ? rkmers2 : kmers2, subref, table);
      if (find_cnt < min_find_cnt)
        continue;
      min_find_cnt = std::max(find_cnt - args.MAX_FIND_CNT_DIFF, min_find_cnt);

      auto sw = SseSmithWaterman::SWResult{};
      if (forward1) {
        if (rprofile2.read == nullptr)
          rprofile2 = SseSmithWaterman::get_profile(rread2);
        sw
          = SseSmithWaterman::align(rprofile2, subref, false, false, min_score);
      } else {
        if (profile2.read == nullptr)
          profile2 = SseSmithWaterman::get_profile(read2);
        sw = SseSmithWaterman::align(profile2, subref, false, false, min_score);
      }

      const auto score = static_cast<int>(sw.score);
      if (score < min_score)
        continue;
      const auto ref_end = sw_pos + sw.ref_end;
      const auto ref_pos = ref_end - sw.read_end;
      SPDLOG_DEBUG("{{ ref pos: {}({}), score: {} }}", ref_pos, static_cast<std::int32_t>(ref_pos-sw_pos), score);

      rescues
        .emplace_back(ref_pos, score, sw.score2, !forward1, sw.read_end,
                      ref_end, find_cnt)
        .rescued
        = true;
      min_score = std::max(min_score, score - args.MAX_SW_DIFF);
    }

    return rescues;
  }

  auto
  shrink_sw_size(int em_size1, auto& sw_alns1, int em_size2,
                 auto& sw_alns2) const {
    if (sw_alns1.size() > args.MAX_SW_CNT)
      sw_alns1.resize(args.MAX_SW_CNT);
    if (em_size1 > args.MAX_EM_CNT)
      sw_alns1.clear();
    if (sw_alns2.size() > args.MAX_SW_CNT)
      sw_alns2.resize(args.MAX_SW_CNT);
    if (em_size2 > args.MAX_EM_CNT)
      sw_alns2.clear();

    if (em_size1 == 0 || em_size2 == 0)
      return;

    const auto sw_size1 = static_cast<int>(sw_alns1.size());
    const auto sw_size2 = static_cast<int>(sw_alns2.size());
    const auto total_size1 = em_size1 + sw_size1;
    const auto total_size2 = em_size2 + sw_size2;

    const auto shrink_size = std::min(total_size1, total_size2);
    if (shrink_size < total_size1) {
      const auto final_sw_size = std::max(0, shrink_size - em_size1);
      if (sw_alns1.size() > final_sw_size)
        sw_alns1.resize(final_sw_size);
    }
    if (shrink_size < total_size2) {
      const auto final_sw_size = std::max(0, shrink_size - em_size2);
      if (sw_alns2.size() > final_sw_size)
        sw_alns2.resize(final_sw_size);
    }
  }

  auto
  insert_penalty(int dist) const {
    const auto ns
      = (dist - args.INSERT_MEAN) / static_cast<double>(args.INSERT_VAR);
    return static_cast<int>(std::pow(ns, 2));
  }

  auto
  get_best_pair(const std::vector<Aln>& alns1, const std::vector<Aln>& alns2,
                const std::vector<AlnPair>& aln_pairs, istring_view read1,
                istring_view rread1, auto& profile1, auto& rprofile1,
                istring_view read2, istring_view rread2, auto& profile2,
                auto& rprofile2, float frac_rep1, float frac_rep2) const {
    SPDLOG_DEBUG("************ pairing results ************");
    print_paires(aln_pairs);

    const auto [opt_score1, sub_score1, sub_cnt1]
      = get_opt_subopt_count(alns1 | std::views::transform(&Aln::score));
    const auto [opt_score2, sub_score2, sub_cnt2]
      = get_opt_subopt_count(alns2 | std::views::transform(&Aln::score));
    const auto [opt_score, sub_score, sub_cnt] = get_opt_subopt_count(
      aln_pairs | std::views::transform(&AlnPair::score));

    auto [aln1, aln2] = aln_pairs.front();
    const auto score_un = opt_score1 + opt_score2 - args.PEN_UNPAIRED;
    const auto success = (opt_score > score_un);
    if (!success) {
      aln1 = alns1.front();
      aln2 = alns2.front();
    }
    if (aln1.forward)
      set_cigar(aln1, read1, profile1);
    else
      set_cigar(aln1, rread1, rprofile1);
    if (aln2.forward)
      set_cigar(aln2, read2, profile2);
    else
      set_cigar(aln2, rread2, rprofile2);
    if (!success) {
      aln1.mapq = mem_approx_mapq_se({aln1.score, aln1.score2, sub_score1,
                                      aln1.align_len, sub_cnt1, frac_rep1});
      aln2.mapq = mem_approx_mapq_se({aln2.score, aln2.score2, sub_score2,
                                      aln2.align_len, sub_cnt2, frac_rep2});
    } else {
      const auto [mapq1, mapq2]
        = mem_mapq_pe({aln1.score, aln1.score2, sub_score1, aln1.align_len,
                       sub_cnt1, frac_rep1},
                      {aln2.score, aln2.score2, sub_score2, aln2.align_len,
                       sub_cnt2, frac_rep2},
                      score_un, opt_score, sub_score, sub_cnt);
      SPDLOG_DEBUG("(raw mapq1:{}, raw mapq2: {})", mapq1, mapq2);
      const auto pen_paired = insert_penalty(aln_pairs.front().dist());
      aln1.mapq = std::max(mapq1 - pen_paired, 0);
      aln2.mapq = std::max(mapq2 - pen_paired, 0);
    }
    aln1.sub_score = (aln1.score == opt_score1 ? sub_score1 : opt_score1);
    aln2.sub_score = (aln2.score == opt_score2 ? sub_score2 : opt_score2);
    return AlnPair{std::move(aln1), std::move(aln2)};
  }

  auto
  map(istring_view read1, istring_view rread1, istring_view read2,
      istring_view rread2) const -> AlnPair {
    SPDLOG_DEBUG("--------------- seeding read1 ---------------");

    const auto [chains1, repeats1, rrepeats1] = seeding(read1, rread1);

    SPDLOG_DEBUG("\n--------------- seeding read2 ---------------");

    const auto [chains2, repeats2, rrepeats2] = seeding(read2, rread2);

    const auto frac_rep1 = (repeats1 + rrepeats1) / (read1.size() * 2.f);
    const auto frac_rep2 = (repeats2 + rrepeats2) / (read2.size() * 2.f);

    auto table = std::vector<bool>(1 << args.KMER_SIZE * 2);
    const auto kmers1 = get_kmers(read1);
    const auto rkmers1 = get_kmers(rread1);
    const auto kmers2 = get_kmers(read2);
    const auto rkmers2 = get_kmers(rread2);

    SPDLOG_DEBUG("\n--------------- exact match read1 ---------------");

    auto [alns1, sw_chains1]
      = exact_match(chains1, read1, rread1, kmers1.size());

    SPDLOG_DEBUG("\n--------------- exact match read2 ---------------");

    auto [alns2, sw_chains2]
      = exact_match(chains2, read2, rread2, kmers2.size());

    SPDLOG_DEBUG("\n--------------- sw read1 ---------------");

    auto [sw_alns1, min_find_cnt1] = get_sw_candidates(
      alns1.empty(), sw_chains1, read1.size(), kmers1, rkmers1, table);

    SPDLOG_DEBUG("\n--------------- sw read2 ---------------");

    auto [sw_alns2, min_find_cnt2] = get_sw_candidates(
      alns2.empty(), sw_chains2, read2.size(), kmers2, rkmers2, table);

    release_memory(sw_chains1);
    release_memory(sw_chains2);

    shrink_sw_size(alns1.size(), sw_alns1, alns2.size(), sw_alns2);

    auto [profile1, rprofile1] = extending(alns1, sw_alns1, read1, rread1);
    auto [profile2, rprofile2] = extending(alns2, sw_alns2, read2, rread2);

    release_memory(sw_alns1);
    release_memory(sw_alns2);

    if (alns1.empty() && alns2.empty()) [[unlikely]]
      return {};

    finalize_alns(alns1);
    finalize_alns(alns2);

    SPDLOG_DEBUG("\n************ force rescue read1 ************");
    auto rescues1 = rescue(alns2, alns1, read1, rread1, profile1, rprofile1,
                           kmers1, rkmers1, table, min_find_cnt1);
    SPDLOG_DEBUG("\n************ force rescue read2 ************");
    auto rescues2 = rescue(alns1, alns2, read2, rread2, profile2, rprofile2,
                           kmers2, rkmers2, table, min_find_cnt2);

    if (!rescues1.empty()) {
      std::ranges::copy(rescues1, std::back_inserter(alns1));
      finalize_alns(alns1);
    }
    SPDLOG_DEBUG("\n************ read1 final result ({}) ************", alns1.size());
    print_alns(alns1);

    if (!rescues2.empty()) {
      std::ranges::copy(rescues2, std::back_inserter(alns2));
      finalize_alns(alns2);
    }
    SPDLOG_DEBUG("\n************ read2 final result ({}) ************", alns2.size());
    print_alns(alns2);

    release_memory(table);
    release_memory(rescues1);
    release_memory(rescues2);

    if (alns2.empty())
      return {
        get_best_one(alns1, read1, rread1, profile1, rprofile1, frac_rep1), {}};
    if (alns1.empty())
      return {
        {}, get_best_one(alns2, read2, rread2, profile2, rprofile2, frac_rep2)};

    SPDLOG_DEBUG("\n--------------- pairing ---------------");

    auto aln_pairs = pairing2(alns1, alns2);
    if (aln_pairs.empty()) {
      SPDLOG_DEBUG("\n--------------- failed ---------------");
      return {
        get_best_one(alns1, read1, rread1, profile1, rprofile1, frac_rep1),
        get_best_one(alns2, read2, rread2, profile2, rprofile2, frac_rep2)};
    }
    return get_best_pair(alns1, alns2, aln_pairs, read1, rread1, profile1,
                         rprofile1, read2, rread2, profile2, rprofile2,
                         frac_rep1, frac_rep2);
  }

 public:
  auto
  map(std::string_view read1, std::string_view read2) const {
    auto iread1 = Codec::to_istring(read1);
    auto iread2 = Codec::to_istring(read2);
    auto riread1 = Codec::rev_comp(iread1);
    auto riread2 = Codec::rev_comp(iread2);
    auto [aln1, aln2] = map(iread1, riread1, iread2, riread2);
    aln1.rev_comp = std::move(riread1);
    aln2.rev_comp = std::move(riread2);
    return std::pair{std::move(aln1), std::move(aln2)};
  }

  auto
  generate_sam(const std::pair<FastqRecord<>, FastqRecord<>>& read) const {
    const auto& name = read.first.name;
    const auto& read1 = read.first.seq;
    const auto& qual1 = read.first.qual;
    const auto& read2 = read.second.seq;
    const auto& qual2 = read.second.qual;

    const auto [aln1, aln2] = map(read1, read2);
    auto [gpos1, score1, score21, forward1, read_end1, ref_end1, find_cnt1,
          align_len1, mapq1, sub_score1, rescued1, cigar1, riread1]
      = aln1;
    auto [gpos2, score2, score22, forward2, read_end2, ref_end2, find_cnt2,
          align_len2, mapq2, sub_score2, rescued2, cigar2, riread2]
      = aln2;
    if (cigar1.empty())
      cigar1 = "*";
    if (cigar2.empty())
      cigar2 = "*";
    // flag1

    auto flag1 = SamUtil::READ_PAIRED + SamUtil::FIRST_OF_PAIR,
         flag2 = SamUtil::READ_PAIRED + SamUtil::SECOND_OF_PAIR;
    if (!forward1) {
      flag1 += SamUtil::READ_REVERSE_STRAND;
      flag2 += SamUtil::MATE_REVERSE_STRAND;
    }
    if (!forward2) {
      flag1 += SamUtil::MATE_REVERSE_STRAND;
      flag2 += SamUtil::READ_REVERSE_STRAND;
    }
    // for WGS:
    // auto [chr1, pos1] = Hs37d5::get_chr_pos(gpos1);
    // auto [chr2, pos2] = Hs37d5::get_chr_pos(gpos2);
    // for single chr
    auto chr = ref.name;
    auto pos1 = gpos1;
    auto pos2 = gpos2;
    auto rname1 = std::string{chr};
    auto rname2 = std::string{chr};
    if (score1 == 0) {
      flag1 += SamUtil::READ_UNMAPPED;
      flag2 += SamUtil::MATE_UNMAPPED;
      rname1 = "*";
    }
    if (score2 == 0) {
      flag1 += SamUtil::MATE_UNMAPPED;
      flag2 += SamUtil::READ_UNMAPPED;
      rname2 = "*";
    }
    auto rnext1 = rname2, rnext2 = rname1;
    auto pnext1 = pos2, pnext2 = pos1;
    auto tlen1 = 0, tlen2 = 0;
    auto proper_pair = false;
    if (score1 != 0 && score2 != 0 && rname1 == rname2) {
      rnext1 = "=";
      rnext2 = "=";
      tlen1
        = SamUtil::compute_tlen(pos1, cigar1, forward1, pos2, cigar2, forward2);
      tlen2 = -tlen1;
      if (forward1 != forward2 && std::abs(tlen1) <= args.PAIR_DIST) {
        flag1 += SamUtil::PROPER_PAIR;
        flag2 += SamUtil::PROPER_PAIR;
        proper_pair = true;
      }
    }

    const auto qname = name.substr(0, name.find_first_of(" \t"));
    auto optionals1 = std::vector<std::string>{
      "AS:i:" + std::to_string(score1), "XS:i:" + std::to_string(sub_score1),
      "RG:Z:1", (rescued1 ? "rs:i:1" : "")};
    auto optionals2 = std::vector<std::string>{
      "AS:i:" + std::to_string(score2), "XS:i:" + std::to_string(sub_score2),
      "RG:Z:1", (rescued2 ? "rs:i:1" : "")};
    auto record1
      = SamRecord{{},  // for base
                  nullptr,
                  qname,
                  static_cast<std::uint16_t>(flag1),
                  rname1,
                  pos1 + 1,
                  mapq1,
                  cigar1,
                  rnext1,
                  pos2 + 1,
                  tlen1,
                  forward1 ? read1 : Codec::to_string(riread1),
                  forward1 ? qual1 : std::string{qual1.rbegin(), qual1.rend()},
                  std::move(optionals1)};

    auto record2
      = SamRecord{{},
                  nullptr,
                  qname,
                  static_cast<std::uint16_t>(flag2),
                  rname2,
                  pos2 + 1,
                  mapq2,
                  cigar2,
                  rnext2,
                  pos1 + 1,
                  tlen2,
                  forward2 ? read2 : Codec::to_string(riread2),
                  forward2 ? qual2 : std::string{qual2.rbegin(), qual2.rend()},
                  std::move(optionals2)};

    return std::pair{std::move(record1), std::move(record2)};
  }
};

}  // namespace biovoltron
