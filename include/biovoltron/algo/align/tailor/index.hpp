#pragma once

#include <biovoltron/algo/align/exact_match/fm_index.hpp>
#include <biovoltron/file_io/fasta.hpp>
#include <biovoltron/utility/istring.hpp>
#include <biovoltron/utility/interval.hpp>
#include <vector>
#include <random>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>

namespace biovoltron {

/**
 * An FM-Index that report chromosome coordinate.
 */
template<
  int SA_INTV = 1,
  typename size_type = std::uint32_t,
  SASorter Sorter = PsaisSorter<size_type>
>
struct Index : FMIndex<SA_INTV,size_type, Sorter> {
  using char_type = std::int8_t;
  
  using Base = FMIndex<SA_INTV, size_type, Sorter>;
  using Base::build;
  using Base::lf;
  using Base::sa_;
  using Base::cnt_;
  using Base::pri_;
  using Base::bwt_;
  using Base::occ_;
  using Base::lookup_;
  using Base::get_offsets;

  Index(int lookup_len = 14) : FMIndex<SA_INTV, size_type, Sorter>{.LOOKUP_LEN = lookup_len} {}
  

  /**
   * Build fm-index by concatening all the chromosome.
   * 
   * @tparam Encoded Whether the sequence is encoded as istring or not.
   * @param ref A list of Fasta file that contains chromosome info.
   * @param lookup_length Length of read for building lookup table.
   */
  template<bool Encoded>
  void make_index(const std::vector<FastaRecord<Encoded>>& ref) {
    static std::random_device r;
    // static std::default_random_engine random_engine(r());
    // std::uniform_int_distribution<ichar> atcg_generator(0, 3);
    chr_bounds.reserve(ref.size());
    auto accu = 0;
    auto ref_seq = istring{};
    for (const auto& record : ref) {
      accu += record.seq.size();
      chr_bounds.emplace_back(record.name, accu-1);
      if constexpr (Encoded)
        ref_seq += record.seq;
      else
        ref_seq += Codec::to_istring(record.seq);
    }
    std::ranges::transform(ref_seq, ref_seq.begin(), [&](auto& c)
      {
        // return c < 4 ? c : atcg_generator(random_engine);
        return c < 4 ? c : 0;
      });
    build(ref_seq);
  }

  /**
   * Get BWT size.
   *
   * @return BWT size.
   */
  auto get_bwt_size() const {
    return bwt_.size();
  }

  /**
   * Transfrom of the index position from the last column to the first column.
   *
   * @param c Character the index position holds.
   * @param i Index of the last column.
   * @return Index of the first column. 
   */
  auto lf_mapping(char_type c, size_type i) const {
    return lf(c, i);
  }
  
  /**
   * Get the size of the input chromosome.
   *
   * @param chr Chromosome name.
   * @throw std::out_of_range if the chromosome is not in the index.
   * @return Chromosome size.
   */
  auto get_chr_size(const std::string_view chr) const {
    auto target_itr = std::ranges::find(chr_bounds, chr, &ChromBound::chrom);
    if (target_itr == chr_bounds.end())
      throw std::out_of_range{"Chromosome is not in the index."};
    if (target_itr == chr_bounds.begin())
      return target_itr->last_elem_pos + 1;
    else
      return target_itr->last_elem_pos - (target_itr-1)->last_elem_pos;
  }

  /**
   * Get exact match results with chromosome coordinate.
   *
   * @param beg Begin of the match range.
   * @param end End of the match range.
   * @param read_len Length of the query read.
   * @return A list of Interval.
   */
  auto get_intervals(size_type beg, size_type end, size_type read_len) const {
    auto intvs = std::vector<Interval>{};
    intvs.reserve(end - beg);

    for (auto pos : get_offsets(beg, end)) {
      auto first = std::ranges::lower_bound(
        chr_bounds, pos, {}, &ChromBound::last_elem_pos);
      auto last = std::ranges::lower_bound(
        chr_bounds, pos+read_len-1, {}, &ChromBound::last_elem_pos);

      // Skip pos out of bound or cross chrom mapping.
      if (first == chr_bounds.end() || first != last)
        continue;

      pos -= (first == chr_bounds.begin()) ? 0 : (first-1)->last_elem_pos + 1;
      intvs.emplace_back(first->chrom, pos, pos + read_len);
    }
    return intvs;
  }

  auto
  save(std::ofstream& fout) const {
    Base::save(fout);

    auto oa = boost::archive::binary_oarchive{fout};
    oa << chr_bounds;
  }

  auto
  load(std::ifstream& fin) {
    Base::load(fin);
    
    auto ia = boost::archive::binary_iarchive{fin};
    ia >> chr_bounds;
    assert(fin.peek() == EOF);
  }

  struct ChromBound {
    std::string chrom;
    std::uint32_t last_elem_pos{};

    friend class boost::serialization::access;

    template<class Archive>
    auto serialize(Archive & ar, const unsigned int version) {
        ar & chrom;
        ar & last_elem_pos;
    }
  };

  std::vector<ChromBound> chr_bounds{};
};

}  // namespace biovoltron
