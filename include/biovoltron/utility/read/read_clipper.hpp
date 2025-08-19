#pragma once
#include <cassert>
#include <biovoltron/file_io/sam.hpp>
#include <biovoltron/utility/interval.hpp>

namespace biovoltron {

/**
 * @ingroup utility
 * @brief Utility class for clipping reads based on CIGAR operations.
 *
 * This class provides methods to hard clip soft-clipped bases, revert
 * soft-clipped bases and hard clip to a specific interval in SAM records.
 */
struct ReadClipper {

  /**
   * Remove bases marked as soft-clipped from both ends of a read.
   *
   * @tparam Encoded Whether the read is encoded (default is false).
   * @param read The SAM record to modify.
   */
  template<bool Encoded = false>
  static void
  hard_clip_soft_clipped_bases(SamRecord<Encoded>& read) {
    auto& seq = read.seq;
    auto& qual = read.qual;
    const auto& cigar = read.cigar;

    auto [front_length, front_op] = cigar.front();
    if (front_op == 'S') {
      seq = seq.substr(front_length);
      qual = qual.substr(front_length);
    }

    auto [back_length, back_op] = cigar.back();
    if (back_op == 'S') {
      seq = seq.substr(0, seq.size() - back_length);
      qual = qual.substr(0, qual.size() - back_length);
    }
  }

  /**
   * Reverts soft-clipped bases at the beginning and end of a read.
   *
   * @tparam Encoded Whether the read is encoded (default is false).
   * @param read The SAM record to modify.
   */
  template<bool Encoded = false>
  static void
  revert_soft_clipped_bases(SamRecord<Encoded>& read) {
    auto& seq = read.seq;
    auto& qual = read.qual;
    auto& cigar = read.cigar;

    if (read.read_reverse_strand()) {
      auto [front_length, front_op] = cigar.front();
      if (front_op == 'S') {
        seq = seq.substr(front_length);
        qual = qual.substr(front_length);
      }
      auto [back_length, back_op] = cigar.back();
      if (back_op == 'S')
        cigar.back() = {back_length, 'M'};
    } else {
      auto [front_length, front_op] = cigar.front();
      auto alignment_begin = read.begin();
      if (front_op == 'S' && alignment_begin >= front_length) {
        cigar.front() = {front_length, 'M'};
        read.pos = alignment_begin - front_length + 1;
      }
      auto [back_length, back_op] = cigar.back();
      if (back_op == 'S') {
        seq = seq.substr(0, seq.size() - back_length);
        qual = qual.substr(0, qual.size() - back_length);
      }
    }
  }

  /**
   * Hard clips a read to fit within a specified interval.
   *
   * @tparam Encoded Whether the read is encoded (default is false).
   * @param read The SAM record to modify.
   * @param interval The interval to hard clip the read to.
   */
  template<bool Encoded = false>
  static void
  hard_clip_to_interval(SamRecord<Encoded>& read, const Interval& interval) {
    auto& seq = read.seq;
    auto& qual = read.qual;

    const auto& [contig, begin, end, strand] = interval;
    assert(read.rname == contig);

    auto alignment_begin = read.begin();
    auto alignment_end = read.end();
    if (alignment_begin < begin) {
      auto clip_size = begin - alignment_begin;
      if (clip_size > seq.size())
        clip_size = seq.size();
      seq = seq.substr(clip_size);
      qual = qual.substr(clip_size);
    }
    if (alignment_end > end) {
      auto clip_size = alignment_end - end;
      seq = seq.substr(0, seq.size() - clip_size);
      qual = qual.substr(0, qual.size() - clip_size);
    }
  }
};

}  // namespace biovoltron
