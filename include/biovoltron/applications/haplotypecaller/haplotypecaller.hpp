#pragma once

#include <biovoltron/algo/align/inexact_match/pairhmm_avx.hpp>
#include <biovoltron/algo/assemble/assembler.hpp>
#include <biovoltron/applications/haplotypecaller/genotype/genotyper.hpp>
#include <biovoltron/file_io/fasta.hpp>
#include <biovoltron/file_io/vcf.hpp>
#include <biovoltron/utility/read/read_clipper.hpp>
#include <biovoltron/utility/read/read_filter.hpp>
#include <random>
#include <spdlog/spdlog.h>
#include <sstream>

namespace biovoltron {

/** 
 * @ingroup applications
 * @brief TBA
 */
struct HaplotypeCaller {
  struct Parameters {
    const int MAX_READS_PER_ALIGN_BEGIN = 5;
    const std::uint32_t REGION_SIZE = 100;
    const std::uint32_t PADDING_SIZE = 100;
  };
  // need transform to upper case
  const FastaRecord<> ref;
  const HaplotypeAssembler assembler;
  const PairHMM pairhmm;
  const Genotyper genotyper;
  const Parameters args;

 private:
  // TODO: this need huge improve
  auto
  generate_reads_map(const std::vector<SamRecord<>>& sam) const {
    std::vector<std::vector<SamRecord<>>> reads_map;
    reads_map.resize(ref.seq.size());
    for (const auto& record : sam)
      if (record.mapq != 0)
        reads_map[record.begin()].emplace_back(record);
    return reads_map;
  }

  auto
  sample_reads(const std::vector<SamRecord<>>& reads) const {
    if (reads.size() <= args.MAX_READS_PER_ALIGN_BEGIN)
      return reads;
    auto sampled_reads = std::vector<SamRecord<>>{};
    std::ranges::sample(reads, std::back_inserter(sampled_reads),
                        args.MAX_READS_PER_ALIGN_BEGIN,
                        std::mt19937{std::random_device{}()});
    return sampled_reads;
  }

  static auto
  filter_reads(std::vector<SamRecord<>>& reads) {
    std::erase_if(reads, MappingQualityReadFilter{});
    std::erase_if(reads, DuplicateReadFilter{});
    std::erase_if(reads, SecondaryAlignmentReadFilter{});
    std::erase_if(reads, MateOnSameContigReadFilter{});
  }

  static auto
  hard_clip_reads(std::vector<SamRecord<>>& reads,
                  const Interval& padded_region) {
    for (auto& read : reads) ReadClipper::revert_soft_clipped_bases(read);
    for (auto& read : reads)
      ReadClipper::hard_clip_to_interval(read, padded_region);
    std::erase_if(reads, MinimumLengthReadFilter{});
  }

  auto
  call_region(std::vector<SamRecord<>>& reads, std::string_view ref,
              const Interval& padded_region,
              const Interval& origin_region) const -> std::vector<Variant> {
    filter_reads(reads);
    hard_clip_reads(reads, padded_region);

    if (reads.empty())
      return {};
    SPDLOG_DEBUG("----------------------------------------------------------------------------------");
    SPDLOG_DEBUG("Assembling {} with {} reads:   (with overlap region = {})", origin_region.to_string(),
      reads.size(), padded_region.to_string());

    auto haplotypes = assembler.assemble(reads, ref);
    if (haplotypes.size() <= 1)
      return {};

    auto likelihoods = pairhmm.compute_likelihoods(haplotypes, reads);
    SPDLOG_DEBUG("----------------------------------------------------------------------------------");
    SPDLOG_DEBUG("Pairhmm values:");
    std::stringstream ss;
    for (auto i = 0; i < likelihoods.size(); i++) {
      ss << "[" << i << "]: ";
      for (auto j = 0; j < likelihoods[i].size(); j++)
        ss << likelihoods[i][j] << " ";
      SPDLOG_DEBUG("{}", ss.str());
      ss.str("");
    }

    SPDLOG_DEBUG("----------------------------------------------------------------------------------");
    SPDLOG_DEBUG("Genotyping:");
    return genotyper.assign_genotype_likelihoods(
      reads, haplotypes, likelihoods, ref, padded_region, origin_region);
  }

 public:
  auto
  run(const std::vector<SamRecord<>>& sam) const {
    const auto ref = static_cast<std::string_view>(this->ref.seq);
    auto raw_variants = std::vector<VcfRecord>{};

    // const auto start_begin = 16050000;
    const auto start_begin = args.PADDING_SIZE;

    auto window_cnt
      = (ref.size() - start_begin + args.REGION_SIZE - 1) / args.REGION_SIZE;
    auto origin_region
      = Interval{this->ref.name, start_begin, start_begin + args.REGION_SIZE};
    auto padded_region = origin_region;
    padded_region.begin -= args.PADDING_SIZE;
    padded_region.end += args.PADDING_SIZE;

    auto reads_map = generate_reads_map(sam);
    for (auto i = 0; i < window_cnt; i++) {
      auto reads = std::vector<SamRecord<>>{};
      for (auto begin = padded_region.begin; begin != padded_region.end;
           begin++)
        if (!reads_map[begin].empty() && begin < reads_map.size()) {
          const auto sampled_reads = sample_reads(reads_map[begin]);
          reads.insert(reads.end(), sampled_reads.begin(), sampled_reads.end());
        }

      if (reads.empty()) {
        SPDLOG_DEBUG("Ignore {}:    (with overlap region = {})", origin_region.to_string(), padded_region.to_string());
      } else {
        const auto region_variants = call_region(
          reads, ref.substr(padded_region.begin, padded_region.size()),
          padded_region, origin_region);
        for (const auto& region_variant : region_variants)
          SPDLOG_DEBUG(region_variant.to_string());
        raw_variants.insert(raw_variants.end(), region_variants.begin(),
                            region_variants.end());
      }

      origin_region.begin += args.REGION_SIZE;
      origin_region.end += args.REGION_SIZE;
      padded_region.begin = origin_region.begin - args.PADDING_SIZE;
      padded_region.end = origin_region.end + args.PADDING_SIZE;
    }
    SPDLOG_DEBUG("HaplotypeCaller done.");
    return raw_variants;
  }
};

}  // namespace biovoltron
