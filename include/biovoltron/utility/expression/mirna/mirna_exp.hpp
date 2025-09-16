#pragma once
#include <biovoltron/file_io/core/record.hpp>
#include <biovoltron/utility/interval.hpp>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <set>
#include <map>
#include <array>

/**
 * @file mirna_exp.hpp
 * @brief Definitions for miRNA expression data structures and operations.
 * 
 * @details This file contains definitions for data structures and operations related to miRNA expression analysis.
 * - @ref biovoltron::mirna::ExpressionMatrix: alias of std::map for expression matrices
 * - @ref biovoltron::mirna::LenExp / TailExp / MirExp: hierarchical expression units
 * - @ref biovoltron::mirna::tail_to_idx and @ref biovoltron::mirna::idx_to_tail: conversion between tail sequence and index
 * - Overloaded arithmetic and aggregation operators
 * 
 * Tail category index conversion: 
 * - 0: 'A', 1: 'C', 2: 'G', 3: 'U', 4: 'O' (Other), 5: 'M' (Genome matching)
 */

 /**
 * @brief Top-level namespace of Biovoltron.
 */
namespace biovoltron {
/**
 * @defgroup mirna_expr miRNA expression data structures
 * @brief Types and utilities for collecting miRNA/iso-miR expression.
 * @{
 */
namespace mirna {
  /**
  * @brief Expression matrix: mapping from miRNA ID to a value type.
  * @tparam T Aggregated expression type (e.g. @ref MirExp).
  *
  * @note Implemented as `std::map<std::string, T>` with ordered keys.
  */
  template<class T>
  using ExpressionMatrix = std::map<std::string, T>;
  
  /**
   * @brief Classify a tail sequence into an index.
   * 
   * @tparam Str String-like type supporting range-based for and substr (e.g. std::string).
   * @param tail Tail sequence string. Empty string denotes genome-matching.
   * 
   * @return category index:
   * - 0: 'A' tail
   * - 1: 'C' tail
   * - 2: 'G' tail
   * - 3: 'U' tail
   * - 4: Other tail (mixed nucleotides)
   * - 5: Genome-matching (empty tail)
   * - -1: Invalid tail (contains non-ACGTU characters)
   * 
   * @note if mutliple nucleotide types are present, returns 4 (Other).
   * @see idx_to_tail
   */
  template <class Str>
  int tail_to_idx(const Str& tail)
  {
    if (tail.empty()) return 5; // genome matching

    std::set<char> nt_check;
    for(auto& nt : tail) nt_check.emplace(nt);

    if(nt_check.size() > 1) return 4; // Other tail

    switch(*nt_check.begin())
    {
      case 'A': return 0; // A tail
      case 'C': return 1; // C tail
      case 'G': return 2; // G tail
      case 'U': case 'T': return 3; // U tail
      default: return -1;
    }
  }

  /**
   * @brief Convert a category index back into a representative character.
   * 
   * @param idx index in range 0..5.
   * @return Representative character:
   * - 0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: 'O' (Other), 5: 'M' (Genome matching)
   * - -1: Invalid index
   * 
   * @see tail_to_idx
   */
  char idx_to_tail(int idx)
  {
    switch(idx)
    {
      case 0: return 'A'; // A tail
      case 1: return 'C'; // C tail
      case 2: return 'G'; // G tail
      case 3: return 'T'; // U tail
      case 4: return 'O'; // Other tail
      case 5: return 'M'; // Genome matching

      default: return -1;
    }
  }
  
  /// @ingroup mirna_expr
  /**
   * @brief Expression contribution for a specific read length.
   * 
   * @details 
   * `LenExp` is the finest-grained unit ; `value` is the expression value for this read length. 
   */
  struct LenExp {
    double value; /// Expression value for this read length.

    /// Add another length-level contribution in place.
    LenExp& operator+=(const LenExp& rhs) {
      value += rhs.value;
      return *this;
    }

    /// Scale the value in place by a constant factor.
    LenExp& operator*=(const double val) {
      value *= val;
      return *this;
    }
  };

  /// @ingroup mirna_expr
  /**
   * @brief Aggregated expression unit for a specific tail category.
   * 
   * @details 
   * - `value`: total expression in this tail category.
   * - `lens`: map from read length to @LenExp
   */
  struct TailExp {
    double value;
    std::map<int, LenExp> lens;
    
    /// Merge another tail-category aggregate.
    TailExp& operator+=(const TailExp& rhs) {
      value += rhs.value;
      for (auto&& [key, value]: rhs.lens) {
          lens[key] += value;
      }
      return *this;
    }

    /// Scale total and all per-length values in place.
    TailExp& operator*=(const double val) {
      value *= val;
      for (auto&& [key, value]: lens) {
          value *= val;
      }
      return *this;
    }
  };
  
  /// @ingroup mirna_expr
  /**
   * @brief Aggregated expression unit for a specific miRNA.
   * 
   * @details
   * - `value`: total expression for this miRNA. (sum of all tails)
   * - `tails`: fixed 6 slots, corresponding to tail categories A/C/G/U/O/M.
   * 
   * @par Hierarchy
   * MirExp -> TailExp[6] -> LenExp map (map<int, LenExp>)
   * 
   * @note Category indices follow @ref tail_to_idx and @ref idx_to_tail
   */
  struct MirExp {
    double value;
    std::array<TailExp, 6> tails;

    /// Aggregate another MirExp (not only value but 6 categories).
    MirExp& operator+=(const MirExp& rhs) {
      value += rhs.value;
      for (auto i = 0; i < tails.size(); ++i) {
        tails[i] += rhs.tails[i];
      }
      return *this;
    }
    /// Scale the expression values by a constant factor.
    MirExp& operator*=(const double val) {
      value *= val;
      for (auto i = 0; i < tails.size(); ++i) {
        tails[i] *= val;
      }
      return *this;
    }

    /**
     * @brief Build a "length * 6 categories" matrix.
     * 
     * @return std::map<int, std::array<double, 6>>
     * - key: read length
     * - value: array of 6 expression values for tail categories A/C/G/U/O/M
     * 
     * @code
     * auto len_based_exp = mir_exp.get_len_based_exp();
     * double a_tail_exp_len_22 = len_based_exp[22][0]; // 22-nt A tail expression
     * @endcode
     */
    std::map<int, std::array<double, 6>> get_len_based_exp() {
      std::map<int, std::array<double, 6>> results;
      for (auto tail_idx = 0; tail_idx < tails.size(); ++tail_idx) {
        const auto& lens_map = tails[tail_idx].lens;
        for (auto&& [len, exp]: lens_map) {
          results[len][tail_idx] += exp.value;
        }
      }
      return results;
    }

    /**
     * @brief Get the partial expression (excluding genome-matching).
     * 
     * @return double Partial expression value (sum of A/C/G/U/O categories).
     * 
     * @note Useful to quantify non-genome-matching miRNA expression.
     */
    double get_partial_exp() {
      double pm = 0.0;
      for (auto tail_idx = 0; tail_idx < 5; ++tail_idx) {
        pm += tails[tail_idx].value;
      }
      return pm;
    }

    /**
     * @brief Initialize MirExp from a single alignment record.
     * 
     * @tparam Aln Alignement-like type providing: 
     * - `aln.seq`: read sequence (string-like)
     * - `aln.hits`: vector of alignment hits (for dilution factor)
     * - `aln.tail_pos`: position of tail start in read sequence (-1 if no tail)
     * @param aln Alignment record.
     * @return MirExp with one category and one length entry.
     * 
     * @details
     * - Dilution factor: `1.0 / aln.hits.size()`
     * - if tail exists: after `aln.tail_pos` is tail, classified by @ref tail_to_idx
     *   and recorded with length `aln.tail_pos`.
     * - if no tail: recorded as genome-matching (index 5) with length `aln.seq.size()`.
     * 
     * @warning Caller must ensure `aln.tail_pos` is valid.
     */
    template <class Aln>
    static MirExp init_from_alignment(Aln&& aln) {
      const auto dilute_factor = 1.0l / aln.hits.size();
      const auto has_tail = aln.tail_pos != -1;
      const double exp = 1.0l * dilute_factor; 
      auto mir_exp = MirExp{};
      auto tail_exp = TailExp{};
      mir_exp.value = exp;

      if (has_tail) {
        const auto read_length = aln.tail_pos; 
        const auto tail_seq = aln.seq.substr(aln.tail_pos);
        const auto tail_idx = tail_to_idx(tail_seq);
        tail_exp.value = exp;
        tail_exp.lens[read_length] = { exp }; 
        mir_exp.tails[tail_idx] = tail_exp;
        return mir_exp;
      } else {
        const auto read_length = aln.seq.size();
        tail_exp.value = exp;
        tail_exp.lens[read_length] = { exp }; 
        // const auto tail_exp = TailExp { exp, {read_length, {exp}} };
        mir_exp.tails[5] = tail_exp;
        return mir_exp;
      }
    } 
  };

  
  /// @relatesalso LenExp
  LenExp operator+(const LenExp& lhs, const LenExp& rhs) {
    return LenExp(lhs) += rhs; 
  }

  /// @relatesalso LenExp
  LenExp operator*(const LenExp& lhs, const double val) {
    return LenExp(lhs) *= val;
  }
  
  /// @relatesalso LenExp
  LenExp operator*(const double val, const LenExp& rhs) {
    return LenExp(rhs) *= val;
  }  
  
  /// @relatesalso TailExp
  TailExp operator+(const TailExp& lhs, const TailExp& rhs) {
    return TailExp(lhs) += rhs; 
  }

  /// @relatesalso TailExp
  TailExp operator*(const TailExp& lhs, const double val) {
    return TailExp(lhs) *= val;
  }
  
  /// @relatesalso TailExp
  TailExp operator*(const double val, const TailExp& rhs) {
    return TailExp(rhs) *= val;
  }  
  
  /// @relatesalso MirExp
  MirExp operator+(const MirExp& lhs, const MirExp& rhs) {
    return MirExp(lhs) += rhs;
  }

  /// @relatesalso MirExp
  MirExp operator*(const MirExp& lhs, const double val) {
    return MirExp(lhs) *= val;
  }
  
  /// @relatesalso MirExp
  MirExp operator*(const double val, const MirExp& rhs) {
    return MirExp(rhs) *= val;
  }
  
/** @} */ // end of group mirna_expr


} // namespace mirna
}  // namespace biovoltron