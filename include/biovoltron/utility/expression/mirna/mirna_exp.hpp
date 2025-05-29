#pragma once
#include <biovoltron/file_io/core/record.hpp>
#include <biovoltron/utility/interval.hpp>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <set>
#include <map>
#include <array>


namespace biovoltron {

namespace mirna {
  template<class T>
  using ExpressionMatrix = std::map<std::string, T>;
  
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

  char idx_to_tail(int idx)
  {
    switch(idx)
    {
      case 0: return 'A'; // A tail
      case 1: return 'C'; // C tail
      case 2: return 'G'; // G tail
      case 3: return 'T'; // U tail
      case 4: return 'O'; // Other tail
      case 5: return 'M'; // Genome mathcing

      default: return -1;
    }
  }
  
  struct LenExp {
    double value;

    LenExp& operator+=(const LenExp& rhs) {
      value += rhs.value;
      return *this;
    }

    LenExp& operator*=(const double val) {
      value *= val;
      return *this;
    }
  };

  struct TailExp {
    double value;
    std::map<int, LenExp> lens;
    
    TailExp& operator+=(const TailExp& rhs) {
      value += rhs.value;
      for (auto&& [key, value]: rhs.lens) {
          lens[key] += value;
      }
      return *this;
    }

    TailExp& operator*=(const double val) {
      value *= val;
      for (auto&& [key, value]: lens) {
          value *= val;
      }
      return *this;
    }
  };
  
  struct MirExp {
    double value;
    std::array<TailExp, 6> tails;

    MirExp& operator+=(const MirExp& rhs) {
      value += rhs.value;
      for (auto i = 0; i < tails.size(); ++i) {
        tails[i] += rhs.tails[i];
      }
      return *this;
    }

    MirExp& operator*=(const double val) {
      value *= val;
      for (auto i = 0; i < tails.size(); ++i) {
        tails[i] *= val;
      }
      return *this;
    }

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

    double get_partial_exp() {
      double pm = 0.0;
      for (auto tail_idx = 0; tail_idx < 5; ++tail_idx) {
        pm += tails[tail_idx].value;
      }
      return pm;
    }

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

  
  
  LenExp operator+(const LenExp& lhs, const LenExp& rhs) {
    return LenExp(lhs) += rhs; 
  }

  LenExp operator*(const LenExp& lhs, const double val) {
    return LenExp(lhs) *= val;
  }
  
  LenExp operator*(const double val, const LenExp& rhs) {
    return LenExp(rhs) *= val;
  }  
  
  TailExp operator+(const TailExp& lhs, const TailExp& rhs) {
    return TailExp(lhs) += rhs; 
  }

  TailExp operator*(const TailExp& lhs, const double val) {
    return TailExp(lhs) *= val;
  }
  
  TailExp operator*(const double val, const TailExp& rhs) {
    return TailExp(rhs) *= val;
  }  
  
  MirExp operator+(const MirExp& lhs, const MirExp& rhs) {
    return MirExp(lhs) += rhs;
  }

  MirExp operator*(const MirExp& lhs, const double val) {
    return MirExp(lhs) *= val;
  }
  
  MirExp operator*(const double val, const MirExp& rhs) {
    return MirExp(rhs) *= val;
  }
  

} // namespace mirna
}  // namespace biovoltron