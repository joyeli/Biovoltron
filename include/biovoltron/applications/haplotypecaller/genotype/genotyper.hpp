#pragma once

#include <biovoltron/applications/haplotypecaller/genotype/allele_frequency/af_calculator.hpp>
#include <biovoltron/applications/haplotypecaller/allele/allele_utils.hpp>
#include <biovoltron/file_io/sam.hpp>
#include <biovoltron/math/math_utils.hpp>
#include <biovoltron/utility/genotype/genotype_utils.hpp>
#include <biovoltron/utility/haplotype/haplotype.hpp>
#include <biovoltron/utility/range/range_utils.hpp>
#include <set>
#include <spdlog/spdlog.h>
#include <sstream>

namespace biovoltron {

class Genotyper {
  static inline const auto SPAN_DEL = std::string{"*"};
  static constexpr auto ALLELE_EXTENSION = 2;
  static constexpr auto MAX_ALLELE_COUNT = GenotypeUtils::MAX_ALLELE_COUNT;
  constexpr static auto STANDARD_CONFIDENCE_FOR_CALLING = 30.0;
  constexpr static auto PHRED_SCALE_QUAL_THRESHOLD
    = STANDARD_CONFIDENCE_FOR_CALLING / -10.0;
  constexpr static auto EPSILON = 1.0E-10;
  constexpr static auto SUM_GL_THRESH_NOCALL = -0.1;

  auto
  process_cigar_for_initial_events(Haplotype& haplotype, std::string_view ref,
                                   const Interval& padded_region) const {
    const auto& cigar = haplotype.cigar;
    const auto& hap = haplotype.seq;
    const auto& [contig, padded_begin, padded_end, strand] = padded_region;

    auto ref_pos = haplotype.align_begin_wrt_ref;
    auto hap_pos = 0;
    for (auto [length, op] : cigar) {
      switch (op) {
        case 'M': {
          auto mismatch_offsets = std::vector<std::uint32_t>{};
          for (auto offset = 0; offset < length; offset++)
            if (ref[ref_pos + offset] != hap[hap_pos + offset])
              mismatch_offsets.push_back(offset);

          if (!mismatch_offsets.empty()) {
            for (auto offset : mismatch_offsets) {
              auto variant = Variant{};
              variant.ref = ref[ref_pos + offset];
              variant.alt = hap[hap_pos + offset];
              const auto mismatch_begin = padded_begin + ref_pos + offset;
              variant.location = {contig, mismatch_begin, mismatch_begin + 1, strand};
              haplotype.event_map.emplace(mismatch_begin, std::move(variant));
            }
          }
          ref_pos += length;
          hap_pos += length;
          break;
        }
        case 'I': {
          if (ref_pos > 0) {
            auto variant = Variant{};
            variant.ref = ref[ref_pos - 1];
            variant.alt = variant.ref + hap.substr(hap_pos, length);
            const auto insertion_begin = padded_begin + ref_pos - 1;
            variant.location = {contig, insertion_begin, insertion_begin + 1, strand};
            haplotype.event_map.emplace(insertion_begin, std::move(variant));
          }
          hap_pos += length;
          break;
        }
        case 'D': {
          if (ref_pos > 0) {
            auto variant = Variant{};
            variant.ref = ref.substr(ref_pos - 1, length + 1);
            variant.alt = ref[ref_pos - 1];
            const auto deletion_begin = padded_begin + ref_pos - 1;
            variant.location
              = {contig, deletion_begin, deletion_begin + length + 1, strand};
            haplotype.event_map.emplace(deletion_begin, std::move(variant));
          }
          ref_pos += length;
          break;
        }
        case 'S': {
          hap_pos += length;
          break;
        }
        default:
          throw std::invalid_argument(
            "Unsupported cigar operator created during SW alignment");
      }
    }
  }

  auto
  set_events_for_haplotypes(std::vector<Haplotype>& haplotypes,
                            std::string_view ref,
                            const Interval& padded_region) const {
    auto events_begins = std::set<std::uint32_t>{};
    auto rank = 0;
    for (auto& h : haplotypes) {
      h.rank = rank++;
      process_cigar_for_initial_events(h, ref, padded_region);
      for (const auto& [begin, even] : h.event_map) events_begins.insert(begin);
    }
    return events_begins;
  }

  auto
  get_events_from_haplotypes(std::uint32_t begin,
                             std::vector<Haplotype>& haplotypes) const {
    auto unique_events = std::set<Variant>{};
    for (const auto& h : haplotypes)
      for (auto& event : h.get_overlapping_events(begin))
        unique_events.insert(std::move(event));
    auto results
      = std::vector<Variant>(unique_events.begin(), unique_events.end());
    return results;
  }

  auto
  replace_span_dels(std::vector<Variant>& events, char ref_allele,
                    std::uint32_t begin) const {
    for (auto& event : events) {
      if (event.location.begin != begin) {
        auto new_event = Variant{};
        new_event.location = {event.location.chrom, begin, begin + 1};
        new_event.ref = ref_allele;
        new_event.alt = SPAN_DEL;
        event = std::move(new_event);
      }
    }
  }

  auto
  determine_reference_allele(const std::vector<Variant>& events) const {
    return std::ranges::max_element(
             events, {}, [](const auto& elem) { return elem.ref.size(); })
      ->ref;
  }

  auto
  get_compatible_alternate_allele(const std::string& ref_allele,
                                  const Variant& event) const {
    if (event.alt == SPAN_DEL)
      return SPAN_DEL;
    return event.alt + ref_allele.substr(event.ref.size());
  }

  auto
  resolve_incompatible_alleles(const std::string& ref_allele,
                               const Variant& event,
                               std::set<std::string>& alts) const {
    if (event.ref == ref_allele)
      alts.insert(event.alt);
    else
      alts.insert(get_compatible_alternate_allele(ref_allele, event));
  }

  auto
  get_compatible_alleles(const std::vector<Variant>& events) const {
    auto longest_event = events.front();
    const auto ref_allele = determine_reference_allele(events);
    auto alleles = std::vector{ref_allele};
    auto alts = std::set<std::string>{};
    for (const auto& event : events) {
      if (event.size() > longest_event.size())
        longest_event = event;
      resolve_incompatible_alleles(ref_allele, event, alts);
    }
    std::ranges::copy(alts, std::back_inserter(alleles));
    return std::pair(alleles, longest_event.location);
  }

  auto
  get_allele_mapper(const std::vector<std::string>& alleles,
                    std::uint32_t begin,
                    const std::vector<Haplotype>& haplotypes) const {
    auto result = std::map<std::uint32_t, std::vector<std::uint32_t>>{};
    result[0];
    const auto& ref_allele = alleles[0];
    for (const auto& h : haplotypes) {
      const auto spanning_events = h.get_overlapping_events(begin);
      if (spanning_events.empty())
        result[0].push_back(h.rank);
      for (const auto& event : spanning_events) {
        if (event.location.begin == begin) {
          if (event.ref.size() == ref_allele.size())
            result[RangeUtils::index_of(alleles, event.alt)].push_back(h.rank);
          else if (event.ref.size() < ref_allele.size())
            result[RangeUtils::index_of(
                     alleles,
                     get_compatible_alternate_allele(ref_allele, event))]
              .push_back(h.rank);
        } else
          result[RangeUtils::index_of(alleles, SPAN_DEL)].push_back(h.rank);
      }
    }
    return result;
  }

  auto
  get_haplotype_mapper(
    const std::map<std::uint32_t, std::vector<std::uint32_t>>& allele_mapper,
    std::uint32_t haplotype_count) const {
    auto haplotype_mapper = std::vector(haplotype_count, 0u);
    for (const auto& [allele_index, haplotype_indices] : allele_mapper)
      for (auto haplotype_index : haplotype_indices)
        haplotype_mapper[haplotype_index] = allele_index;
    return haplotype_mapper;
  }

  auto
  get_read_indices_to_keep(const std::vector<SamRecord<>>& reads,
                           const Interval& overlap) const {
    auto read_indices_to_keep = std::vector<std::uint32_t>{};
    read_indices_to_keep.reserve(reads.size());
    for (auto i = 0; i < reads.size(); i++)
      if (static_cast<Interval>(reads[i]).overlaps(overlap))
        read_indices_to_keep.push_back(i);
    return read_indices_to_keep;
  }

  auto
  marginal_likelihoods(
    std::uint32_t allele_count,
    const std::vector<std::uint32_t>& haplotype_mapper,
    const std::vector<std::uint32_t>& read_indices_to_keep,
    const std::vector<std::vector<double>>& haplotype_likelihoods) const {
    auto allele_likelihoods = std::vector(
      read_indices_to_keep.size(),
      std::vector(allele_count, std::numeric_limits<double>::lowest()));
    for (auto r = 0; r < read_indices_to_keep.size(); r++) {
      const auto old_read_index = read_indices_to_keep[r];
      for (auto h = 0; h < haplotype_mapper.size(); h++) {
        const auto allele_index = haplotype_mapper[h];
        const auto likelihood = haplotype_likelihoods[old_read_index][h];
        if (likelihood > allele_likelihoods[r][allele_index])
          allele_likelihoods[r][allele_index] = likelihood;
      }
    }
    return allele_likelihoods;
  }

  auto
  marginalize(const std::vector<std::uint32_t>& haplotype_mapper,
              std::uint32_t allele_count, const std::vector<SamRecord<>>& reads,
              const std::vector<std::vector<double>>& haplotype_likelihoods,
              const Interval& overlap) const {
    const auto read_indices_to_keep = get_read_indices_to_keep(reads, overlap);
    return marginal_likelihoods(allele_count, haplotype_mapper,
                                read_indices_to_keep, haplotype_likelihoods);
  }

  auto
  single_component_genotype_likelihood_by_read(
    std::vector<double>& genotype_likelihoods,
    const std::vector<std::vector<double>>& allele_likelihoods,
    std::uint32_t a) const {
    const auto log10_frequency = std::log10(2);
    std::ranges::transform(allele_likelihoods,
                           std::back_inserter(genotype_likelihoods),
                           [=](const auto& likelihoods) {
                             return likelihoods[a] + log10_frequency;
                           });
  }

  auto
  two_component_genotype_likelihood_by_read(
    std::vector<double>& genotype_likelihoods,
    const std::vector<std::vector<double>>& allele_likelihoods,
    std::uint32_t a1, std::uint32_t a2) const {
    std::ranges::transform(allele_likelihoods,
                           std::back_inserter(genotype_likelihoods),
                           [=](const auto& likelihoods) {
                             return MathUtils::approximate_log10_sum_log10(
                               likelihoods[a1], likelihoods[a2]);
                           });
  }

  auto
  calculate_read_likelihoods_by_genotype_index(
    const std::vector<std::vector<double>>& allele_likelihoods,
    std::uint32_t allele_count) const {
    auto read_likelihoods_by_genotype_index = std::vector(
      (allele_count + 1) * allele_count / 2, std::vector<double>{});
    auto cur_genotype_index = 0;
    for (auto a1 = 0; a1 < allele_count; a1++) {
      for (auto a2 = a1; a2 < allele_count; a2++) {
        auto& read_genotype_likelihoods
          = read_likelihoods_by_genotype_index[cur_genotype_index++];
        read_genotype_likelihoods.reserve(allele_likelihoods.size());
        if (a1 == a2)
          single_component_genotype_likelihood_by_read(
            read_genotype_likelihoods, allele_likelihoods, a1);
        else
          two_component_genotype_likelihood_by_read(read_genotype_likelihoods,
                                                    allele_likelihoods, a1, a2);
      }
    }
    return read_likelihoods_by_genotype_index;
  }

  auto
  get_genotype_likelihoods(const std::vector<std::vector<double>>&
                             read_likelihoods_by_genotype_index) const {
    const auto genotype_count = read_likelihoods_by_genotype_index.size();
    auto result = std::vector(genotype_count, 0.0);
    const auto denominator
      = read_likelihoods_by_genotype_index[0].size() * std::log10(2);
    for (auto genotype = 0; genotype < genotype_count; genotype++)
      result[genotype]
        = std::accumulate(read_likelihoods_by_genotype_index[genotype].begin(),
                          read_likelihoods_by_genotype_index[genotype].end(),
                          0.0)
          - denominator;
    return result;
  }

  auto
  calculate_genotype_likelihoods(
    const std::vector<std::vector<double>>& allele_likelihoods,
    std::uint32_t allele_count) const {
    const auto read_likelihoods_by_genotype_index
      = calculate_read_likelihoods_by_genotype_index(allele_likelihoods,
                                                     allele_count);
    return get_genotype_likelihoods(read_likelihoods_by_genotype_index);
  }

  static auto
  calculate_output_allele_subset(
    std::span<const std::string> alleles,
    const std::vector<std::pair<std::string, std::pair<int, double>>>&
      mle_counts_and_log10p_ref_by_allele) {
    auto output_alleles = std::vector<std::pair<std::string, int>>{};
    auto site_is_monomorphic = true;
    for (const auto& allele : alleles | std::views::drop(1)) {
      const auto& [mle_count, log10p_non_ref]
        = std::ranges::find(
            mle_counts_and_log10p_ref_by_allele | std::views::keys, allele)
            .base()
            ->second;
      const auto is_plausible
        = (log10p_non_ref + EPSILON) < PHRED_SCALE_QUAL_THRESHOLD;

      site_is_monomorphic &= !is_plausible;
      if (is_plausible)
        output_alleles.emplace_back(allele, mle_count);
    }
    return std::pair{std::move(output_alleles), site_is_monomorphic};
  }

  static auto
  calculate_genotypes(std::span<const double> log10_genotype_likelihoods,
                      std::span<const std::string> alleles)
    -> std::tuple<std::vector<std::string>, double, std::vector<double>> {
    const auto genotypes = GenotypeUtils::get_vcf_genotypes(alleles.size());
    const auto [log10p_no_variant, mle_counts_and_log10p_ref_by_allele]
      = AFCalculator::calculate(log10_genotype_likelihoods, alleles, genotypes);
    const auto [output_alleles, site_is_monomorphic]
      = calculate_output_allele_subset(alleles,
                                       mle_counts_and_log10p_ref_by_allele);
    if (output_alleles.empty())
      return {};

    // note the math.abs is necessary because -10 * 0.0 => -0.0 which isn't nice
    const auto log10_confidence
      = !site_is_monomorphic ?
          log10p_no_variant + 0.0 :
          MathUtils::log10_one_minus_pow10(log10p_no_variant) + 0.0;

    const auto phred_scaled_confidence = (-10.0 * log10_confidence) + 0.0;

    if (phred_scaled_confidence < STANDARD_CONFIDENCE_FOR_CALLING)
      return {};
    if (output_alleles.size() == 1 && output_alleles.begin()->first == SPAN_DEL)
      return {};

    auto resulting_alleles = std::vector{alleles[0]};
    std::ranges::copy(output_alleles | std::views::keys,
                      std::back_inserter(resulting_alleles));

    const auto new_likelihoods = AlleleUtils::subset_alleles(
      log10_genotype_likelihoods, alleles, resulting_alleles, genotypes);
    if (std::accumulate(new_likelihoods.begin(), new_likelihoods.end(), 0.0)
        >= SUM_GL_THRESH_NOCALL)
      return {};

    return {std::move(resulting_alleles), phred_scaled_confidence,
            std::move(new_likelihoods)};
  }

 public:
  auto
  assign_genotype_likelihoods(
    const std::vector<SamRecord<>>& reads, std::vector<Haplotype>& haplotypes,
    const std::vector<std::vector<double>>& haplotype_likelihoods,
    std::string_view ref, const Interval& padded_region,
    const Interval& origin_region) const {
    const auto events_begins
      = set_events_for_haplotypes(haplotypes, ref, padded_region);

    SPDLOG_DEBUG("events_begins:");
    for (const auto& begin : events_begins) SPDLOG_DEBUG("{}", begin);

    const auto& [contig, origin_begin, origin_end, strand] = origin_region;
    auto variants = std::vector<Variant>{};
    for (auto begin : events_begins) {
      if (begin < origin_begin || begin >= origin_end)
        continue;
      auto events = get_events_from_haplotypes(begin, haplotypes);

      SPDLOG_DEBUG("events:");
      for (const auto& event : events) SPDLOG_DEBUG("{}", event.to_string());

      replace_span_dels(events, ref[begin - padded_region.begin], begin);
      auto [alleles, alleles_loc] = get_compatible_alleles(events);

      SPDLOG_DEBUG("alleles:");
      std::stringstream ss;
      for (const auto& allele : alleles) ss << allele << " ";
      SPDLOG_DEBUG("{}", ss.str());

      const auto allele_count = alleles.size();
      if (allele_count > MAX_ALLELE_COUNT)
        continue;
      const auto allele_mapper = get_allele_mapper(alleles, begin, haplotypes);
      const auto haplotype_mapper
        = get_haplotype_mapper(allele_mapper, haplotypes.size());
      const auto allele_likelihoods = marginalize(
        haplotype_mapper, allele_count, reads, haplotype_likelihoods,
        alleles_loc.expand_with(ALLELE_EXTENSION));
      const auto genotype_likelihoods
        = calculate_genotype_likelihoods(allele_likelihoods, allele_count);

      if (std::ranges::max_element(genotype_likelihoods)
          == genotype_likelihoods.begin())
        continue;

      const auto pls = GenotypeUtils::gls_to_pls(genotype_likelihoods);
      const auto log10_genotype_likelihoods = [&pls] {
        auto new_gls = std::vector<double>{};
        for (const auto pl : pls) new_gls.push_back(pl / -10.0);
        return new_gls;
      }();

      SPDLOG_DEBUG("genotype_likelihoods:");
      ss.str("");
      for (const auto& likelihoods : genotype_likelihoods)
        ss << likelihoods << " ";
      SPDLOG_DEBUG("{}", ss.str());

      SPDLOG_DEBUG("pls:");
      ss.str("");
      for (const auto& pl : pls) ss << pl << " ";
      SPDLOG_DEBUG("{}", ss.str());

      SPDLOG_DEBUG("log10_genotype_likelihoods:");
      ss.str("");
      for (const auto& likelihoods : log10_genotype_likelihoods)
        ss << likelihoods << " ";
      SPDLOG_DEBUG("{}", ss.str());

      const auto vcf_genotype_likelihoods
        = GenotypeUtils::to_vcf_order(log10_genotype_likelihoods);

      const auto [output_alleles, phred_scaled_confidence, new_likelihoods]
        = calculate_genotypes(vcf_genotype_likelihoods, alleles);

      SPDLOG_DEBUG("phred_scaled_confidence: {}", phred_scaled_confidence);

      SPDLOG_DEBUG("output_alleles:");
      ss.str("");
      for (const auto& allele : output_alleles) ss << allele << " ";
      SPDLOG_DEBUG("{}", ss.str());

      SPDLOG_DEBUG("new_likelihoods:");
      ss.str("");
      for (const auto& likelihoods : new_likelihoods)
        ss << likelihoods << " ";
      SPDLOG_DEBUG("{}", ss.str());

      if (output_alleles.empty())
        continue;

      const auto genotype_index
        = std::ranges::max_element(new_likelihoods) - new_likelihoods.begin();
      const auto genotype = GenotypeUtils::get_vcf_genotypes(
        output_alleles.size())[genotype_index];

      const auto new_pls = [&new_likelihoods] {
        auto new_pls = std::vector<int>{};
        for (const auto likelihoods : new_likelihoods)
          new_pls.push_back(likelihoods * -10.0);
        return new_pls;
      }();

      variants.push_back({.location = alleles_loc,
                          .alleles = output_alleles,
                          .gt = genotype,
                          .pls = new_pls,
                          .gq = RangeUtils::second<std::ranges::less>(new_pls),
                          .qual = phred_scaled_confidence});
      SPDLOG_DEBUG("{}", variants.back().to_string());
    }
    return variants;
  }
};

}  // namespace biovoltron
