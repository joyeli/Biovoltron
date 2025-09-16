// Copyright (c) 2020 Robert Vaser

#ifndef SIMD_ALIGNMENT_ENGINE_HPP_
#define SIMD_ALIGNMENT_ENGINE_HPP_

#include <cstdint>
#include <memory>
#include <string>
#include <vector>
#include <exception>
#include <limits>
#include <stdexcept>
#include <iostream>
#include <biovoltron/algo/align/spoa/graph.hpp>
#include <simdpp/simd.h>




namespace biovoltron{

  /**
    * This project is revised from Spoa (https://github.com/rvaser/spoa) 
    * Spoa (SIMD POA) is a C++ implementation of the partial order alignment (POA) algorithm
    * which is used to generate consensus sequences. It supports three alignment modes: local
    * (Smith-Waterman), global (Needleman-Wunsch), and semi-global alignment (overlap), and three gap
    * modes: linear and affine. It also supports Intel SSE4.1+ and AVX2 vectorization.
    *
    * Example
    * ```cpp
    * #include <iostream>
    * #include "Biovoltron/include/biovoltron/algo/align/spoa/simd_alignment_engine.hpp"
    * #include "Biovoltron/include/biovoltron/algo/align/spoa/graph.hpp"
    * 
    * int main(int argc, char** argv) {
    * 
    *   std::vector<std::string> sequences = {
    *       "CATAAAAGAACGTAGGTCGCCCGTCCGTAACCTGTCGGATCACCGGAAAGGACCCGTAAAGTGATAATGAT",
    *       "ATAAAGGCAGTCGCTCTGTAAGCTGTCGATTCACCGGAAAGATGGCGTTACCACGTAAAGTGATAATGATTAT",
    *       "ATCAAAGAACGTGTAGCCTGTCCGTAATCTAGCGCATTTCACACGAGACCCGCGTAATGGG",
    *       "CGTAAATAGGTAATGATTATCATTACATATCACAACTAGGGCCGTATTAATCATGATATCATCA",
    *       "GTCGCTAGAGGCATCGTGAGTCGCTTCCGTACCGCAAGGATGACGAGTCACTTAAAGTGATAAT",
    *       "CCGTAACCTTCATCGGATCACCGGAAAGGACCCGTAAATAGACCTGATTATCATCTACAT"
    *   };
    * 
    *   auto alignment_engine = spoa::AlignmentEngine::Create(
    *       spoa::AlignmentType::kNW, 3, -5, -3);  // linear gaps
    * 
    *   spoa::Graph graph{};
    * 
    *   for (const auto& it : sequences) {
    *     auto alignment = alignment_engine->Align(it, graph);
    *     graph.AddAlignment(alignment, it);
    *   }
    * 
    *   auto consensus = graph.GenerateConsensus();
    * 
    *   std::cerr << ">Consensus LN:i:" << consensus.size() << std::endl
    *             << consensus << std::endl;
    * 
    *   auto msa = graph.GenerateMultipleSequenceAlignment();
    * 
    *   for (const auto& it : msa) {
    *     std::cerr << it << std::endl;
    *   }
    * 
    *   return 0;
    * }
    * ```
    */

#ifdef __AVX2__  
  /**
    * T_NUM times paralelism of simd
    */
  constexpr std::size_t T_NUM = 16;

  constexpr std::size_t LSS = 2;
  constexpr std::size_t RSS = 30;

  /**
    * In prefix max, it check the previous cell in the same row
    * Similar to "H_row[j] = std::max(H_row[j - 1] + g_, H_row[j])" in single thread.
    */
  template <int N>
  inline void prefix_max(simdpp::int16<16>& a,
                        const simdpp::int16<16> (&penalties)[N],
                        const simdpp::int16<16> (&masks)[N]) {
      simdpp::int16<16> m;
      m = simdpp::add(a, penalties[0]);
      a = simdpp::max(a, masks[0] | simdpp::move8_r<1>(m));

      m = simdpp::add(a, penalties[1]);
      a = simdpp::max(a, masks[1] | simdpp::move8_r<2>(m));

      m = simdpp::add(a, penalties[2]);
      a = simdpp::max(a, masks[2] | simdpp::move8_r<4>(m));

      m = simdpp::add(a, penalties[3]);
      a = simdpp::max(a, masks[3] | simdpp::move8_r<8>(m));
  }

#elif defined(__SSE2__)
  /**
    * Specifies the level of parallelism for SIMD operations.
    */
  constexpr std::size_t T_NUM = 8;


  constexpr std::size_t LSS = 2;
  constexpr std::size_t RSS = 14;

  /**
   * In prefix max, it compare the cell to its left cell.
   * Similar to "H_row[j] = std::max(H_row[j - 1] + g_, H_row[j])" in single thread.
   */
  template <int N>
  inline void prefix_max(simdpp::int16<8>& a,
                        const simdpp::int16<8> (&penalties)[N],
                        const simdpp::int16<8> (&masks)[N]) {
      simdpp::int16<8> m;
      m = simdpp::add(a, penalties[0]);
      a = simdpp::max(a, masks[0] | simdpp::move8_r<1>(m));

      m = simdpp::add(a, penalties[1]);
      a = simdpp::max(a, masks[1] | simdpp::move8_r<2>(m));

      m = simdpp::add(a, penalties[2]);
      a = simdpp::max(a, masks[2] | simdpp::move8_r<4>(m));
  }
#else
  #error "Unsupported SIMD instruction set"
#endif

constexpr std::int16_t kNegativeInfinity =
    std::numeric_limits<std::int16_t>::min() + 1024;


enum class AlignmentType {
  kSW,  // Smith Waterman
  kNW   // Needleman Wunsch
};

enum class AlignmentSubtype {
  kLinear,  // g * i
  kAffine
};

class Graph;
using Alignment = std::vector<std::pair<std::int32_t, std::int32_t>>;

class SimdAlignmentEngine{
 private:
    AlignmentType type_;
    AlignmentSubtype subtype_;
    std::int8_t m_;
    std::int8_t n_;
    std::int8_t g_;
    std::int8_t e_;
    std::int8_t q_;
    std::int8_t c_;

    struct Implementation {
        std::vector<std::uint32_t> node_id_to_rank;
        std::vector<std::int32_t> M;
        std::int16_t* sequence_profile;
        std::int16_t* H;
        std::int16_t* F;
        std::int16_t* E;
        std::int32_t* O;
        std::int32_t* Q;
        std::int32_t last_length;

        Implementation()
            : node_id_to_rank(),
                sequence_profile(),
                M(),
                H(nullptr),
                F(nullptr),
                E(nullptr),
                O(nullptr),
                Q(nullptr),
                last_length(0) {
        }
    };

    std::unique_ptr<Implementation> pimpl_;


    SimdAlignmentEngine(
    AlignmentType type,
    AlignmentSubtype subtype,
    std::int8_t m,
    std::int8_t n,
    std::int8_t g,
    std::int8_t e,
    std::int8_t q,
    std::int8_t c)
    : type_(type),
      subtype_(subtype),
      m_(m),
      n_(n),
      g_(g),
      e_(e),
      q_(q),
      c_(c),
      pimpl_(new SimdAlignmentEngine::Implementation()) {}


 public:
  SimdAlignmentEngine(const SimdAlignmentEngine&) = delete;
  SimdAlignmentEngine& operator=(const SimdAlignmentEngine&) = delete;

  SimdAlignmentEngine(SimdAlignmentEngine&&) = default;
  SimdAlignmentEngine& operator=(SimdAlignmentEngine&&) = default;

  ~SimdAlignmentEngine() = default;

  static std::unique_ptr<SimdAlignmentEngine> Create(
    AlignmentType type,
    std::int8_t m,   // match
    std::int8_t n,   // mismatch
    std::int8_t g){
      return Create(type, m, n, g, g);
      };

  static std::unique_ptr<SimdAlignmentEngine> Create(
      AlignmentType type,
      std::int8_t m,
      std::int8_t n,
      std::int8_t g,   // gap open
      std::int8_t e    // gap extention
      ){
        return Create(type, m, n, g, e, g, e);
        };
    
  static std::unique_ptr<SimdAlignmentEngine> Create(
      AlignmentType type,
      std::int8_t m,
      std::int8_t n,
      std::int8_t g,
      std::int8_t e,
      std::int8_t q,
      std::int8_t c) {
    if (type != AlignmentType::kSW &&
        type != AlignmentType::kNW) {
        throw std::invalid_argument(
            "[spoa::AlignmentEngine::Create] error: invalid alignment type!");
    }
    if (g > 0 || q > 0) {
        throw std::invalid_argument(
            "[spoa::AlignmentEngine::Create] error: "
            "gap opening penalty must be non-positive!");
    }
    if (e > 0 || c > 0) {
        throw std::invalid_argument(
            "[spoa::AlignmentEngine::Create] error: "
            "gap extension penalty must be non-positive!");
    }

    AlignmentSubtype subtype = g >= e ?
        AlignmentSubtype::kLinear :  AlignmentSubtype::kAffine;

    if (subtype == AlignmentSubtype::kLinear) {
        e = g;
    } else if (subtype == AlignmentSubtype::kAffine) {
        q = g;
        c = e;
    }
    return std::unique_ptr<SimdAlignmentEngine>(
      new SimdAlignmentEngine(type, subtype, m, n, g, e, q, c));
    
  };

  static std::unique_ptr<SimdAlignmentEngine> Create(
      AlignmentType type,
      AlignmentSubtype subtype,
      std::int8_t m,
      std::int8_t n,
      std::int8_t g,
      std::int8_t e,
      std::int8_t q,
      std::int8_t c);

  /**
    * Check if the resource is enough. If the resource is enough, call Reallocate to allocate the space of DP table.
    */
  void Prealloc(
      std::uint32_t max_sequence_len,
      std::uint8_t alphabet_size) {
    if (max_sequence_len > std::numeric_limits<int32_t>::max()) {
      throw std::invalid_argument(
          "[spoa::SimdAlignmentEngine::Prealloc] error: too large sequence!");
    }
    try {
      Realloc(
          static_cast<std::uint64_t>(max_sequence_len) + 1,
          static_cast<std::uint64_t>(max_sequence_len) * alphabet_size + alphabet_size,  // NOLINT
          alphabet_size);
    } catch (std::bad_alloc& ba) {
      throw std::invalid_argument(
          "[spoa::SimdAlignmentEngine::Prealloc] error: insufficient memory!");
    }
  }


  /**
    * allocate the DP table
    */
  void Realloc(
      std::uint64_t matrix_width,
      std::uint64_t matrix_height,
      std::uint8_t num_codes) {

    std::uint64_t allc_width = (matrix_width / T_NUM + 1) * T_NUM;

    if (pimpl_->node_id_to_rank.size() < matrix_height - 1) {
      pimpl_->node_id_to_rank.resize(matrix_height - 1, 0);
    }


    if (subtype_ == AlignmentSubtype::kLinear) {
      if (pimpl_->last_length < matrix_height * allc_width) {
        pimpl_->last_length = matrix_height * allc_width * 3;
        pimpl_->sequence_profile = (std::int16_t*) aligned_alloc(32, pimpl_->last_length * sizeof(std::int16_t));
        pimpl_->H = (std::int16_t*) aligned_alloc(32, pimpl_->last_length * sizeof(std::int16_t));
        pimpl_->F = nullptr;
        pimpl_->E = nullptr;
      }
    }  
    else if (subtype_ == AlignmentSubtype::kAffine) {
      if (pimpl_->last_length < matrix_height * allc_width) {
        pimpl_->last_length = matrix_height * allc_width * 9;
        pimpl_->sequence_profile = (std::int16_t*) aligned_alloc(32, pimpl_->last_length * sizeof(std::int16_t));
        pimpl_->H = (std::int16_t*) aligned_alloc(32, pimpl_->last_length * sizeof(std::int16_t));
        pimpl_->F = (std::int16_t*) aligned_alloc(32, pimpl_->last_length * sizeof(std::int16_t));
        pimpl_->E = (std::int16_t*) aligned_alloc(32, pimpl_->last_length * sizeof(std::int16_t));
      }
    }
  }

  void Initialize(
      const char* sequence, std::uint32_t sequence_len,
      const Graph& graph) noexcept {
    std::uint32_t matrix_width = sequence_len + 1;
    std::uint32_t matrix_height = graph.nodes().size() + 1;
    std::uint32_t allc_width = (matrix_width / T_NUM + 1) * T_NUM;

    // initialize the table of match score
    for (std::uint32_t i = 0; i < graph.num_codes(); ++i) {
      char c = graph.decoder(i);
      pimpl_->sequence_profile[i * allc_width] = 0;
      for (std::uint32_t j = 0; j < sequence_len; ++j) {
        pimpl_->sequence_profile[i * allc_width + (j + 1)] =
            (c == sequence[j] ? m_ : n_);
      }
    }

    const auto& rank_to_node = graph.rank_to_node();
    for (std::uint32_t i = 0; i < rank_to_node.size(); ++i) {
      pimpl_->node_id_to_rank[rank_to_node[i]->id] = i;
    }

    // initialize secondary matrices
    switch (subtype_) {
      case AlignmentSubtype::kLinear:
        pimpl_->H[0] = 0;
        break;
      case AlignmentSubtype::kAffine:
        pimpl_->F[0] = 0;
        pimpl_->E[0] = 0;
        for (std::uint32_t j = 1; j < matrix_width; ++j) {
          pimpl_->F[j] = kNegativeInfinity;
          pimpl_->E[j] = g_ + (j - 1) * e_;
        }
        for (std::uint32_t i = 1; i < matrix_height; ++i) {
          const auto& edges = rank_to_node[i - 1]->inedges;
          std::int16_t penalty = edges.empty() ? g_ - e_ : kNegativeInfinity;
          for (const auto& it : edges) {
            std::uint32_t pred_i = pimpl_->node_id_to_rank[it->tail->id] + 1;
            penalty = std::max(penalty, pimpl_->F[pred_i * allc_width]);
          }
          pimpl_->F[i * allc_width] = penalty + e_;
          pimpl_->E[i * allc_width] = kNegativeInfinity;
        }
      default:
        break;
    }

    // initialize primary matrix
    switch (type_) {
      case AlignmentType::kSW:
        for (std::uint32_t j = 1; j < allc_width; ++j) {
          pimpl_->H[j] = 0;
        }
        for (std::uint32_t i = 1; i < matrix_height; ++i) {
          pimpl_->H[i * allc_width] = 0;
        }
        break;
      case AlignmentType::kNW:
        switch (subtype_) {
          case AlignmentSubtype::kLinear:
            for (std::uint32_t j = 1; j < allc_width; ++j) {
              pimpl_->H[j] = j * g_;
            }
            for (std::uint32_t i = 1; i < matrix_height; ++i) {
              const auto& edges = rank_to_node[i - 1]->inedges;
              std::int16_t penalty = edges.empty() ? 0 : kNegativeInfinity;
              for (const auto& it : edges) {
                std::uint32_t pred_i = pimpl_->node_id_to_rank[it->tail->id] + 1;
                penalty = std::max(penalty, pimpl_->H[pred_i * allc_width]);
              }
              pimpl_->H[i * allc_width] = penalty + g_;
            }
            break;
          case AlignmentSubtype::kAffine:
            for (std::uint32_t j = 1; j < matrix_width; ++j) {
              pimpl_->H[j] = pimpl_->E[j];
            }
            for (std::uint32_t i = 1; i < matrix_height; ++i) {
              pimpl_->H[i * allc_width] = pimpl_->F[i * allc_width];
            }
            break;
          default:
            break;
        }
        break;
      default:
        break;
    }
  }

  /**
    * Alignment engine dispatcher
    */
  Alignment Align(
      const char* sequence, std::uint32_t sequence_len,
      const Graph& graph,
      std::int32_t* score) {

      if (sequence_len > std::numeric_limits<int32_t>::max()) {
          throw std::invalid_argument(
              "[spoa::SimdAlignmentEngine::Align] error: too large sequence!");
      }

      if (graph.nodes().empty() || sequence_len == 0) {
          return Alignment();
      }

      if (WorstCaseAlignmentScore(sequence_len, graph.nodes().size()) < kNegativeInfinity) {  // NOLINT
          throw std::invalid_argument(
              "[spoa::SimdAlignmentEngine::Align] error: possible overflow!");
      }

      try {
          Realloc(sequence_len + 1, graph.nodes().size() + 1, graph.num_codes());
      } catch (std::bad_alloc& ba) {
          throw std::invalid_argument(
              "[spoa::SimdAlignmentEngine::Align] error: insufficient memory!");
      }
      Initialize(sequence, sequence_len, graph);

      if (subtype_ == AlignmentSubtype::kLinear) {
          return Linear(sequence_len, graph, score);
      } else if (subtype_ == AlignmentSubtype::kAffine) {
          return Affine(sequence_len, graph, score);
      } 
      return Alignment();
  }

  Alignment Align(
      const std::string& sequence,
      const Graph& graph,
      std::int32_t* score = nullptr) {
  return Align(sequence.c_str(), sequence.size(), graph, score);
  }


  std::int64_t WorstCaseAlignmentScore(
      std::int64_t i,
      std::int64_t j) const {
  auto gap_score = [&] (std::int64_t len) -> std::int64_t {
      return len == 0 ? 0 : std::min(g_ + (len - 1) * e_, q_ + (len - 1) * c_);
  };
  return std::min(
      -1 * (m_ * std::min(i, j) + gap_score(std::abs(i - j))),
      gap_score(i) + gap_score(j));
  }


private:
  /**
    * Linear gap alignment
    */
  Alignment Linear(
      std::uint32_t sequence_len,
      const Graph& graph,
      std::int32_t* score) noexcept {
    std::uint64_t matrix_width = sequence_len + 1;
    std::uint64_t allc_width = (matrix_width / T_NUM + 1) * T_NUM;
    const auto& rank_to_node = graph.rank_to_node();

    std::int16_t max_score = type_ == AlignmentType::kSW ? 0 : kNegativeInfinity;
    std::uint32_t max_i = 0;
    std::uint32_t max_j = 0;
    auto update_max_score = [&max_score, &max_i, &max_j] (
        std::int16_t* H_row,
        std::uint32_t i,
        std::uint32_t j) -> void {
      if (max_score < H_row[j]) {
        max_score = H_row[j];
        max_i = i;
        max_j = j;
      }
      return;
    };

    // initialization
    simdpp::int16<T_NUM> penalties[5];
    simdpp::int16<T_NUM> masks[5];

    // pentalties
    penalties[0] = simdpp::make_int(g_);
    for (std::uint32_t i = 1; i < 4; ++i) {
      penalties[i] = simdpp::add(penalties[i - 1], penalties[i - 1]);
    }

    // masks
    __attribute__((aligned(32))) int16_t unpacked[16] = {0};  
    for (std::uint32_t i = 0, j = 0; i < T_NUM; ++i) {
      unpacked[i] = kNegativeInfinity;
      if ((i & (i + 1)) == 0) {
        masks[j++] = simdpp::load(unpacked);
      }
    }

    masks[4] = simdpp::load(unpacked);
    masks[4] = simdpp::move8_r<LSS/2>(masks[4]);
    

    // first
    simdpp::int16<T_NUM> gap = simdpp::make_int(g_);
    simdpp::int16<T_NUM> kNegSimd = simdpp::make_int(kNegativeInfinity);
    simdpp::int16<T_NUM> x, tmp, a;

    /**
      * The representation of POA (Partial Order Alignment) graph view.
      *
      * ```
      *       A  C  G  C  A  A
      *  A    .  .  .  .  .  .
      *  | \
      *  A T
      *  | /
      *  C
      *  |
      *  C
      *  |
      *  A
      *  | \
      *  | G
      *  |/
      *  A
      *  ```
      *
      * Takes the final node A for example, the predecessor of A includes of A & G.
      * Each row represents a node in the POA graph, and each column corresponds to the alignment
      * of a new read to the POA graph. 
      */
    

    // alignment
    for (const auto& it : rank_to_node) {
      const auto& char_profile =
          &(pimpl_->sequence_profile[it->code * allc_width]);

      std::uint32_t i = pimpl_->node_id_to_rank[it->id] + 1;
      std::uint32_t pred_i = it->inedges.empty() ?
          0 : pimpl_->node_id_to_rank[it->inedges[0]->tail->id] + 1;

      std::int16_t* H_row = &(pimpl_->H[i * allc_width]);
      std::int16_t* H_pred_row = &(pimpl_->H[pred_i * allc_width]);

      // update H
      for (std::uint64_t j = 1; j < matrix_width; ++j) {
        H_row[j] = std::max(
            H_pred_row[j - 1] + char_profile[j],
            H_pred_row[j] + g_);
      }

      // check other predeccessors
      for (std::uint32_t p = 1; p < it->inedges.size(); ++p) {
        pred_i = pimpl_->node_id_to_rank[it->inedges[p]->tail->id] + 1;

        H_pred_row = &(pimpl_->H[pred_i * allc_width]);

        for (std::uint64_t j = 1; j < matrix_width; ++j) {
          H_row[j] = std::max(
              (int16_t) (H_pred_row[j - 1] + char_profile[j]),
              std::max(
                  H_row[j],
                  (int16_t) (H_pred_row[j] + g_)));
        }
      }

      x = simdpp::move8_l<RSS/2>(kNegSimd);
      for (std::uint64_t j = 0; j < matrix_width; j += T_NUM) {
        // Check the previous cell in the same row, implementing "H_row[j] = std::max(static_cast<int16_t>(H_row[j - 1] + g), H_row[j])" as in a scalar (SISD) context.
        a = simdpp::load(H_row + j);
        a = simdpp::max(a, x | masks[4]);

        prefix_max(a, penalties, masks);

        tmp = simdpp::add(a, gap);
        x = tmp >> LSS;
        simdpp::store(H_row + j, a);
      }

      if(it->outedges.empty()){
        update_max_score(H_row, i, matrix_width-1);
      }


    }

    if (max_i == 0 && max_j == 0) {
      return Alignment();
    }
    if (score) {
      *score = max_score;
    }

    // backtrack
    Alignment alignment;
    std::uint32_t i = max_i;
    std::uint32_t j = max_j;

    auto sw_condition = [this, &i, &j, &allc_width] () -> bool {
      return (pimpl_->H[i * allc_width + j] == 0) ? false : true;
    };
    auto nw_condition = [&i, &j] () -> bool {
      return (i == 0 && j == 0) ? false : true;
    };
    auto ov_condition = [&i, &j] () -> bool {
      return (i == 0 || j == 0) ? false : true;
    };

    std::uint32_t prev_i = 0;
    std::uint32_t prev_j = 0;
    while ((type_ == AlignmentType::kSW && sw_condition()) ||
          (type_ == AlignmentType::kNW && nw_condition())) {
      auto H_ij = pimpl_->H[i * allc_width + j];
      bool predecessor_found = false;

      if (i != 0 && j != 0) {
        const auto& it = rank_to_node[i - 1];
        std::int32_t match_cost =
            pimpl_->sequence_profile[it->code * allc_width + j];

        std::uint32_t pred_i = it->inedges.empty() ?
            0 : pimpl_->node_id_to_rank[it->inedges[0]->tail->id] + 1;

        if (H_ij == pimpl_->H[pred_i * allc_width + (j - 1)] + match_cost) {
          prev_i = pred_i;
          prev_j = j - 1;
          predecessor_found = true;
        } else {
          for (std::uint32_t p = 1; p < it->inedges.size(); ++p) {
            std::uint32_t pred_i =
                pimpl_->node_id_to_rank[it->inedges[p]->tail->id] + 1;

            if (H_ij == pimpl_->H[pred_i * allc_width + (j - 1)] + match_cost) {
              prev_i = pred_i;
              prev_j = j - 1;
              predecessor_found = true;
              break;
            }
          }
        }
      }

      if (!predecessor_found && i != 0) {
        const auto& it = rank_to_node[i - 1];

        std::uint32_t pred_i = it->inedges.empty() ? 0 :
            pimpl_->node_id_to_rank[it->inedges[0]->tail->id] + 1;

        if (H_ij == pimpl_->H[pred_i * allc_width + j] + g_) {
          prev_i = pred_i;
          prev_j = j;
          predecessor_found = true;
        } else {
          for (std::uint32_t p = 1; p < it->inedges.size(); ++p) {
            std::uint32_t pred_i =
                pimpl_->node_id_to_rank[it->inedges[p]->tail->id] + 1;

            if (H_ij == pimpl_->H[pred_i * allc_width + j] + g_) {
              prev_i = pred_i;
              prev_j = j;
              predecessor_found = true;
              break;
            }
          }
        }
      }

      if (!predecessor_found && H_ij == pimpl_->H[i * allc_width + j - 1] + g_) {  // NOLINT
        prev_i = i;
        prev_j = j - 1;
        predecessor_found = true;
      }

      alignment.emplace_back(
          i == prev_i ? -1 : rank_to_node[i - 1]->id,
          j == prev_j ? -1 : j - 1);

      i = prev_i;
      j = prev_j;
    }

    std::reverse(alignment.begin(), alignment.end());
    return alignment;
  }




  Alignment Affine(
      std::uint32_t sequence_len,
      const Graph& graph,
      std::int32_t* score) noexcept {

    std::uint64_t matrix_width = sequence_len + 1;
    std::uint64_t allc_width = (matrix_width / T_NUM + 1) * T_NUM;
    const auto& rank_to_node = graph.rank_to_node();

    std::int16_t max_score = type_ == AlignmentType::kSW ? 0 : kNegativeInfinity;
    std::uint32_t max_i = 0;
    std::uint32_t max_j = 0;
    auto update_max_score = [&max_score, &max_i, &max_j] (
        std::int16_t* H_row,
        std::uint32_t i,
        std::uint32_t j) -> void {
      if (max_score < H_row[j]) {
        max_score = H_row[j];
        max_i = i;
        max_j = j;
      }
      return;
    };

    // initialization
    simdpp::int16<T_NUM> penalties[5];
    simdpp::int16<T_NUM> masks[5];

    // pentalties
    penalties[0] = simdpp::make_int(e_);
    for (std::uint32_t i = 1; i < 4; ++i) {
      penalties[i] = simdpp::add(penalties[i - 1], penalties[i - 1]);
    }

    // masks
    __attribute__((aligned(32))) int16_t unpacked[T_NUM] = {0};  
    for (std::uint32_t i = 0, j = 0; i < T_NUM; ++i) {
      unpacked[i] = kNegativeInfinity;
      if ((i & (i + 1)) == 0) {
        masks[j++] = simdpp::load(unpacked);
      }
    }

    // add
    simdpp::int16<T_NUM> H_row_simd, F_row_simd, E_row_simd;
    simdpp::int16<T_NUM> extend = simdpp::make_int(e_);
    simdpp::int16<T_NUM> gap = simdpp::make_int(g_ - e_);
    simdpp::int16<T_NUM> x, m, a;

    // alignment
    for (const auto& it : rank_to_node) {
      const auto& char_profile =
          &(pimpl_->sequence_profile[it->code * allc_width]);

      std::uint32_t i = pimpl_->node_id_to_rank[it->id] + 1;
      std::uint32_t pred_i = it->inedges.empty() ? 0 :
          pimpl_->node_id_to_rank[it->inedges[0]->tail->id] + 1;

      std::int16_t* H_row = &(pimpl_->H[i * allc_width]);
      std::int16_t* H_pred_row = &(pimpl_->H[pred_i * allc_width]);

      std::int16_t* F_row = &(pimpl_->F[i * allc_width]);
      std::int16_t* F_pred_row = &(pimpl_->F[pred_i * allc_width]);

      // update F and H
      for (std::uint64_t j = 1; j < matrix_width; ++j) {
        F_row[j] = std::max(
            H_pred_row[j] + g_,
            F_pred_row[j] + e_);
        H_row[j] = H_pred_row[j - 1] + char_profile[j];
      }
      // check other predeccessors
      for (std::uint32_t p = 1; p < it->inedges.size(); ++p) {
        pred_i = pimpl_->node_id_to_rank[it->inedges[p]->tail->id] + 1;

        H_pred_row = &(pimpl_->H[pred_i * allc_width]);
        F_pred_row = &(pimpl_->F[pred_i * allc_width]);

        for (std::uint64_t j = 1; j < matrix_width; ++j) {
          F_row[j] = std::max(
              F_row[j],
              (int16_t) std::max(
                  H_pred_row[j] + g_,
                  F_pred_row[j] + e_));
          H_row[j] = std::max(
              H_row[j],
              (int16_t) (H_pred_row[j - 1] + char_profile[j]));
        }
      }

      // update E and H
      std::int16_t* E_row = &(pimpl_->E[i * allc_width]);
      x = simdpp::make_int(kNegativeInfinity);
      for (std::uint64_t j = 0; j < matrix_width; j += T_NUM) {
        H_row_simd = simdpp::load(H_row + j);
        F_row_simd = simdpp::load(F_row + j);

        H_row_simd = simdpp::max(H_row_simd, F_row_simd);

        E_row_simd = simdpp::add(
          simdpp::add(
            simdpp::move8_r<LSS/2>(H_row_simd) | simdpp::move8_l<RSS/2>(x),
            gap
          ),
          extend
        );

        prefix_max(E_row_simd, penalties, masks);
        H_row_simd = simdpp::max(H_row_simd, E_row_simd);

        x = simdpp::max(H_row_simd, simdpp::sub(E_row_simd, gap));
        
        simdpp::store(E_row + j, E_row_simd);
        simdpp::store(H_row + j, H_row_simd);

        if(it->outedges.empty()){
          update_max_score(H_row, i, matrix_width-1);
        }
      }
    }

    if (max_i == 0 && max_j == 0) {
      return Alignment();
    }
    if (score) {
      *score = max_score;
    }

    // backtrack
  Alignment alignment;
  std::uint32_t i = max_i;
  std::uint32_t j = max_j;

  auto sw_condition = [this, &i, &j, &allc_width] () -> bool {
    return (pimpl_->H[i * allc_width + j] == 0) ? false : true;
  };
  auto nw_condition = [&i, &j] () -> bool {
    return (i == 0 && j == 0) ? false : true;
  };
  auto ov_condition = [&i, &j] () -> bool {
    return (i == 0 || j == 0) ? false : true;
  };

  std::uint32_t prev_i = 0;
  std::uint32_t prev_j = 0;

  while ((type_ == AlignmentType::kSW && sw_condition()) ||
         (type_ == AlignmentType::kNW && nw_condition())) {
    auto H_ij = pimpl_->H[i * allc_width + j];
    bool predecessor_found = false, extend_left = false, extend_up = false;

    if (i != 0 && j != 0) {
      const auto& it = rank_to_node[i - 1];
      std::int32_t match_cost =
          pimpl_->sequence_profile[it->code * allc_width + j];

      std::uint32_t pred_i = it->inedges.empty() ?
          0 : pimpl_->node_id_to_rank[it->inedges[0]->tail->id] + 1;

      if (H_ij == pimpl_->H[pred_i * allc_width + (j - 1)] + match_cost) {
        prev_i = pred_i;
        prev_j = j - 1;
        predecessor_found = true;
      } else {
        for (std::uint32_t p = 1; p < it->inedges.size(); ++p) {
          pred_i = pimpl_->node_id_to_rank[it->inedges[p]->tail->id] + 1;

          if (H_ij == pimpl_->H[pred_i * allc_width + (j - 1)] + match_cost) {
            prev_i = pred_i;
            prev_j = j - 1;
            predecessor_found = true;
            break;
          }
        }
      }
    }

    if (!predecessor_found && i != 0) {
      const auto& it = rank_to_node[i - 1];

      std::uint32_t pred_i = it->inedges.empty() ?
          0 : pimpl_->node_id_to_rank[it->inedges[0]->tail->id] + 1;

      if ((extend_up = H_ij == pimpl_->F[pred_i * allc_width + j] + e_) ||
                       H_ij == pimpl_->H[pred_i * allc_width + j] + g_) {
        prev_i = pred_i;
        prev_j = j;
        predecessor_found = true;
      } else {
        for (std::uint32_t p = 1; p < it->inedges.size(); ++p) {
          pred_i = pimpl_->node_id_to_rank[it->inedges[p]->tail->id] + 1;

          if ((extend_up = H_ij == pimpl_->F[pred_i * allc_width + j] + e_) ||
                           H_ij == pimpl_->H[pred_i * allc_width + j] + g_) {
            prev_i = pred_i;
            prev_j = j;
            predecessor_found = true;
            break;
          }
        }
      }
    }

    if (!predecessor_found && j != 0) {
      if ((extend_left = H_ij == pimpl_->E[i * allc_width + j - 1] + e_) ||
                         H_ij == pimpl_->H[i * allc_width + j - 1] + g_) {
        prev_i = i;
        prev_j = j - 1;
        predecessor_found = true;
      }
    }

    alignment.emplace_back(
        i == prev_i ? -1 : rank_to_node[i - 1]->id,
        j == prev_j ? -1 : j - 1);

    i = prev_i;
    j = prev_j;

    if (extend_left) {
      while (true) {
        alignment.emplace_back(-1, j - 1);
        --j;
        if (pimpl_->E[i * allc_width + j] + e_ !=
            pimpl_->E[i * allc_width + j + 1]) {
          break;
        }
      }
    } else if (extend_up) {
      while (true) {
        bool stop = false;
        prev_i = 0;
        for (const auto& it : rank_to_node[i - 1]->inedges) {
          std::uint32_t pred_i = pimpl_->node_id_to_rank[it->tail->id] + 1;

          if ((stop = pimpl_->F[i * allc_width + j] == pimpl_->H[pred_i * allc_width + j] + g_) ||  // NOLINT
                      pimpl_->F[i * allc_width + j] == pimpl_->F[pred_i * allc_width + j] + e_) {  // NOLINT
            prev_i = pred_i;
            break;
          }
        }

        alignment.emplace_back(rank_to_node[i - 1]->id, -1);
        i = prev_i;
        if (stop || i == 0) {
          break;
        }
      }
    }
  }

  std::reverse(alignment.begin(), alignment.end());
    return alignment;
  }


};

}  // namespace spoa

#endif  // SIMD_ALIGNMENT_ENGINE_HPP_
