#ifndef WFPOA_ALIGNMENT_ENGINE_HPP_
#define WFPOA_ALIGNMENT_ENGINE_HPP_

// #define WFUNIT

#include <cstdint>
#include <memory>
#include <string>
#include <vector>
#include <exception>
#include <limits>
#include <stdexcept>
#include <string>
#include <iostream>
#include <algorithm>
#include <unordered_set>
#include <queue>
#include <chrono> 
#include <biovoltron/algo/align/wfpoa/graph.hpp>


namespace biovoltron{

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


class WfGraph;
using Alignment = std::vector<std::pair<std::int32_t, std::int32_t>>;

class WfpoaAlignmentEngine{
 private:
    AlignmentType type_;
    AlignmentSubtype subtype_;
    std::int8_t m_;
    std::int8_t n_;
    std::int8_t g_;
    std::int8_t e_;
    std::int8_t q_;
    std::int8_t c_;

    int d, cut_threshold = std::numeric_limits<int>::max(), record;
    int mxi, mxj;
    std::uint32_t matrix_width;
    std::string seq;
    std::string code_to_char{""};
    std::queue<std::pair<int, int>> S;


    struct Implementation {
      std::vector<std::uint32_t> node_id_to_rank;
      std::vector<std::int32_t> sequence_profile;
      std::vector<std::int32_t> M;
      std::vector<std::vector<int>> i_to_next_i;
      std::int32_t* H;
      std::int32_t* F;
      std::int32_t* E;
      std::int32_t* O;
      std::int32_t* Q;

      Implementation()
          : node_id_to_rank(),
            sequence_profile(),
            M(),
            i_to_next_i(),
            H(nullptr),
            F(nullptr),
            E(nullptr),
            O(nullptr),
            Q(nullptr) {
      }
    };

    std::unique_ptr<Implementation> pimpl_;


    WfpoaAlignmentEngine(
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
      pimpl_(new WfpoaAlignmentEngine::Implementation()) {

      total_time_measure_spoa_fulltable = 0;
      total_time_measure_spoa_wavefront = 0;
      }


 public:
    bool acc_bench = false;
    double total_time_measure_spoa_fulltable = 0;
    double total_time_measure_spoa_wavefront = 0;
  WfpoaAlignmentEngine(const WfpoaAlignmentEngine&) = delete;
  WfpoaAlignmentEngine& operator=(const WfpoaAlignmentEngine&) = delete;

  WfpoaAlignmentEngine(WfpoaAlignmentEngine&&) = default;
  WfpoaAlignmentEngine& operator=(WfpoaAlignmentEngine&&) = default;

  ~WfpoaAlignmentEngine() = default;

  void setCutThreshold(int cut){
    cut_threshold = cut;
  }

  static std::unique_ptr<WfpoaAlignmentEngine> Create(
    AlignmentType type,
    std::int8_t m,   // match
    std::int8_t n,   // mismatch
    std::int8_t g){
      return Create(type, m, n, g, g);
      };

  static std::unique_ptr<WfpoaAlignmentEngine> Create(
      AlignmentType type,
      std::int8_t m,
      std::int8_t n,
      std::int8_t g,   // gap open
      std::int8_t e    // gap extention
      ){
        return Create(type, m, n, g, e, g, e);
        };
    
    static std::unique_ptr<WfpoaAlignmentEngine> Create(
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
      return std::unique_ptr<WfpoaAlignmentEngine>(
        new WfpoaAlignmentEngine(type, subtype, m, n, g, e, q, c));
      
    };

    static std::unique_ptr<WfpoaAlignmentEngine> Create(
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
            "[spoa::WfpoaAlignmentEngine::Prealloc] error: too large sequence!");
      }
      try {
        Realloc(
            static_cast<std::uint64_t>(max_sequence_len) + 1,
            static_cast<std::uint64_t>(max_sequence_len) * alphabet_size + alphabet_size,  // NOLINT
            alphabet_size);
      } catch (std::bad_alloc& ba) {
        throw std::invalid_argument(
            "[spoa::WfpoaAlignmentEngine::Prealloc] error: insufficient memory!");
      }
    }


    /**
      * allocate the DP table
      */
    void Realloc(
        std::uint64_t matrix_width,
        std::uint64_t matrix_height,
        std::uint8_t num_codes) {
      if (pimpl_->node_id_to_rank.size() < matrix_height - 1) {
        pimpl_->node_id_to_rank.resize(matrix_height - 1, 0);
      }
      if (pimpl_->sequence_profile.size() < num_codes * matrix_width) {
        pimpl_->sequence_profile.resize(num_codes * matrix_width, 0);
      }
      if (subtype_ == AlignmentSubtype::kLinear) {
        if (pimpl_->M.size() < matrix_height * matrix_width) {
          pimpl_->M.resize(matrix_width * matrix_height, 0);
          pimpl_->H = pimpl_->M.data();
          pimpl_->F = nullptr;
          pimpl_->E = nullptr;
        }
      } else if (subtype_ == AlignmentSubtype::kAffine) {
        if (pimpl_->M.size() < 3 * matrix_height * matrix_width) {
          pimpl_->M.resize(3 * matrix_width * matrix_height, 0);
          pimpl_->H = pimpl_->M.data();
          pimpl_->F = pimpl_->H + matrix_width * matrix_height;
          pimpl_->E = pimpl_->F + matrix_width * matrix_height;
        }
      } 
    }

  void wf_Realloc_Init(const WfGraph& graph, int matrix_height){      
      const auto& rank_to_node = graph.rank_to_node();

      for (std::uint32_t i = 0; i < graph.num_codes(); ++i) {
        code_to_char.push_back(graph.decoder(i));
      }

      // const auto& rank_to_node = graph.rank_to_node();
      for (std::uint32_t i = 0; i < rank_to_node.size(); ++i) {
        pimpl_->node_id_to_rank[rank_to_node[i]->id] = i;
      }

      pimpl_->i_to_next_i.resize(matrix_height, {});
      pimpl_->i_to_next_i[0].clear();
      for(auto &n : graph.first_node()){
        pimpl_->i_to_next_i[0].push_back(n+1);
      }

      int i = 1;
      for(auto& cur_node : rank_to_node){
        pimpl_->i_to_next_i[i].clear();
        for(auto& edge : cur_node->outedges){
          pimpl_->i_to_next_i[i].push_back(pimpl_->node_id_to_rank[edge->head->id]+1);
        }
        i += 1;
      }
  }


  void Initialize(
      const char* sequence, std::uint32_t sequence_len,
      const WfGraph& graph) noexcept {
    std::uint32_t matrix_width = sequence_len + 1;
    std::uint32_t matrix_height = graph.nodes().size() + 1;

    for (std::uint32_t i = 0; i < graph.num_codes(); ++i) {
      char c = graph.decoder(i);
      code_to_char.push_back(c);
      pimpl_->sequence_profile[i * matrix_width] = 0;
      for (std::uint32_t j = 0; j < sequence_len; ++j) {
        pimpl_->sequence_profile[i * matrix_width + (j + 1)] =
            (c == sequence[j] ? m_ : n_);
      }
    }

    const auto& rank_to_node = graph.rank_to_node();
    for (std::uint32_t i = 0; i < rank_to_node.size(); ++i) {
      pimpl_->node_id_to_rank[rank_to_node[i]->id] = i;
    }

    // initialize secondary matrices
    switch (subtype_) {
      case AlignmentSubtype::kAffine:
        pimpl_->F[0] = 0;
        pimpl_->E[0] = 0;
        for (std::uint32_t j = 1; j < matrix_width; ++j) {
          pimpl_->F[j] = kNegativeInfinity;
          pimpl_->E[j] = g_ + (j - 1) * e_;
        }
        for (std::uint32_t i = 1; i < matrix_height; ++i) {
          const auto& edges = rank_to_node[i - 1]->inedges;
          std::int32_t penalty = edges.empty() ? g_ - e_ : kNegativeInfinity;
          for (const auto& it : edges) {
            std::uint32_t pred_i = pimpl_->node_id_to_rank[it->tail->id] + 1;
            penalty = std::max(penalty, pimpl_->F[pred_i * matrix_width]);
          }
          pimpl_->F[i * matrix_width] = penalty + e_;
          pimpl_->E[i * matrix_width] = kNegativeInfinity;
        }
        // fall through
      case AlignmentSubtype::kLinear:
        pimpl_->H[0] = 0;
        break;
      default:
        break;
    }

    // initialize primary matrix
    switch (type_) {
      case AlignmentType::kSW:
        for (std::uint32_t j = 1; j < matrix_width; ++j) {
          pimpl_->H[j] = 0;
        }
        for (std::uint32_t i = 1; i < matrix_height; ++i) {
          pimpl_->H[i * matrix_width] = 0;
        }
        break;
      case AlignmentType::kNW:
        switch (subtype_) {
          case AlignmentSubtype::kAffine:
            for (std::uint32_t j = 1; j < matrix_width; ++j) {
              pimpl_->H[j] = pimpl_->E[j];
            }
            for (std::uint32_t i = 1; i < matrix_height; ++i) {
              pimpl_->H[i * matrix_width] = pimpl_->F[i * matrix_width];
            }
            break;
          case AlignmentSubtype::kLinear:
            for (std::uint32_t j = 1; j < matrix_width; ++j) {
              pimpl_->H[j] = j * g_;
            }
            for (std::uint32_t i = 1; i < matrix_height; ++i) {
              const auto& edges = rank_to_node[i - 1]->inedges;
              std::int32_t penalty = edges.empty() ? 0 : kNegativeInfinity;
              for (const auto& it : edges) {
                std::uint32_t pred_i = pimpl_->node_id_to_rank[it->tail->id] + 1;
                penalty = std::max(penalty, pimpl_->H[pred_i * matrix_width]);
              }
              pimpl_->H[i * matrix_width] = penalty + g_;
            }
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
      const WfGraph& graph,
      std::int32_t* score, 
      int method) {

      if (sequence_len > std::numeric_limits<int32_t>::max()) {
          throw std::invalid_argument(
              "[spoa::WfpoaAlignmentEngine::Align] error: too large sequence!");
      }

      if (graph.nodes().empty() || sequence_len == 0) {
          return Alignment();
      }

      if (WorstCaseAlignmentScore(sequence_len, graph.nodes().size()) < kNegativeInfinity) {  // NOLINT
          throw std::invalid_argument(
              "[spoa::WfpoaAlignmentEngine::Align] error: possible overflow!");
      }



      try {
          Realloc(sequence_len + 1, graph.nodes().size() + 1, graph.num_codes());
      } catch (std::bad_alloc& ba) {
          throw std::invalid_argument(
              "[spoa::WfpoaAlignmentEngine::Align] error: insufficient memory!");
      }
      // Initialize(sequence, sequence_len, graph);

      wf_Realloc_Init(graph, graph.nodes().size() + 1);

// #ifdef showcompare
//       return Debug(sequence_len, graph, score, sequence);
// #endif

      return WFAlignment(sequence_len, graph, score, sequence);
      // return Debug(sequence_len, graph, score, sequence);

      // if(method == 1)
      //   return Linear(sequence_len, graph, score, sequence);

      // return WFAlignment(sequence_len, graph, score, sequence);
  
  }

  Alignment Align(
      const std::string& sequence,
      const WfGraph& graph,
      int method) {
  return Align(sequence.c_str(), sequence.size(), graph, nullptr, method);
  }



  Alignment Align(
      const std::string& sequence,
      const WfGraph& graph,
      int method, std::int32_t* score) {
  return Align(sequence.c_str(), sequence.size(), graph, score, method);
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

    void print_M(int *M, std::vector<WfGraph::Node*> rank_to_node){

      std::cout << "\t\t";
      for(auto c : seq){
          std::cout << "\t\t" << c;
      }
      std::cout << std::endl;

      for (int i = 0; i <= rank_to_node.size(); i ++){
        if(i == 0){
          std::cout << "\t\t";
        }else{
          std::cout << code_to_char[rank_to_node[i-1]->code] << "\t\t";
        }

        for (int j = 0; j < matrix_width; j ++){
          std::cout << M[i * matrix_width + j] << "\t\t";
        }
        std::cout << std::endl;
      }

      std::cout << std::endl << std::endl;

      std::cout << code_to_char[(*(rank_to_node.end() - 1))->code] << ' ' << (*(rank_to_node.end() - 1))->id << std::endl;



    }

    bool extend_arrow(std::vector<int> &M, const WfGraph& graph){
      const auto& rank_to_node = graph.rank_to_node();
      const auto& i_to_next_i = pimpl_->i_to_next_i;
      std::queue<std::pair<int, int>> S_;

      while(!S.empty()){
          auto [i, j] = S.front();
          S.pop();

          record = std::max(record, j);

          if(j < record - cut_threshold) continue;

          if(j == matrix_width -1 ){
            mxi = i;
            mxj = j;
            return false;
          }

          if(i_to_next_i[i].size() == 0){
            if(j < matrix_width - 1){
              S_.push({i, j});
              continue;
            }
          }

          for(auto& next_i : i_to_next_i[i]){
            int pos = next_i * matrix_width + j;
            if(code_to_char[rank_to_node[next_i-1]->code] != seq[j]){
              S_.push({i, j});
            }else if(M[pos + 1] > 0 && M[pos] > 0 && M[pos + 1 - matrix_width] > 0){
              M[pos + 1] = d;
              S.push({next_i, j + 1});
            }
          }

      }


      S = std::move(S_);
      return true;
    }

    bool extend(std::vector<int> &M, const WfGraph& graph){
      const auto& rank_to_node = graph.rank_to_node();
      const auto& i_to_next_i = pimpl_->i_to_next_i;
      std::queue<std::pair<int, int>> S_;

      while(!S.empty()){
          auto [i, j] = S.front();
          S.pop();
          S_.push({i, j});

          record = std::max(record, j);

          if(j < record - cut_threshold) {
            continue;
          }

          if(i_to_next_i[i].size() == 0 && j == matrix_width - 1){
              mxi = i;
              mxj = j;
            return false;
          }

          for(auto &next_i : i_to_next_i[i]){
              if(code_to_char[rank_to_node[next_i-1]->code] == seq[j] && M[next_i * matrix_width + j + 1] > 0){
                  M[next_i * matrix_width + j + 1] = d;
                  S.push({next_i, j+1});
              }
          }
      }


      S = std::move(S_);
      return true;
    }

    void expand(std::vector<int> &M, const WfGraph& graph){
      const auto& rank_to_node = graph.rank_to_node();
      const auto& i_to_next_i = pimpl_->i_to_next_i;
      std::queue<std::pair<int, int>> S_;

        while(!S.empty()){
            auto [i, j] = S.front();
            S.pop();

            if(j == matrix_width || j < record - cut_threshold){
              continue;
            }

            bool is_edge = (j == matrix_width -1);

            if(!is_edge && M[i * matrix_width + j + 1]  > 0){
              M[i * matrix_width + j + 1] = d;
              S_.push({i, j + 1});
            }

            for(auto& next_i : i_to_next_i[i]){
                auto pos = next_i * matrix_width + j;
                if(M[pos] > 0){
                    S_.push({next_i, j});
                    M[pos] = d;
                }

#ifdef WFUNIT
                if(!is_edge && M[pos+1] > 0){
                    S_.push({next_i, j + 1});
                    M[pos+1] = d;
                }
#endif
            }
        }

      S = std::move(S_);
    }

    Alignment WFAlignment(
        std::uint32_t sequence_len,
        const WfGraph& graph,
        std::int32_t* score,
        const char* sequence) noexcept {
      std::uint64_t matrix_width = sequence_len + 1;
      const auto& rank_to_node = graph.rank_to_node();

      std::int32_t max_score = type_ == AlignmentType::kSW ? 0 : kNegativeInfinity;
      std::uint32_t max_i = 0;
      std::uint32_t max_j = 0;


     std::vector<int> M(matrix_width * (rank_to_node.size() + 1), 257);

      std::string seq(sequence);
      this->seq = seq;    
      this->matrix_width = matrix_width;

      S = {};
      S.push({0, 0});
      M[0] = 0;
      d = record = 0;


      auto start_ = std::chrono::high_resolution_clock::now();


      while(
#ifdef WFARROW
        extend_arrow(M, graph)
#else
        extend(M, graph)
#endif
        ){
          d--;     
          expand(M, graph);
      } 

      auto end_ = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double, std::micro> duration_ = end_ - start_;
      total_time_measure_spoa_wavefront += static_cast<double>(duration_.count() );


      if (mxi == 0 && mxj == 0) {
        return Alignment();
      }
      if (score) {
        *score = max_score;
      }

      auto al = traceback(mxi, mxj, rank_to_node, matrix_width, M.data());

      return al;
    }


    Alignment Linear(
        std::uint32_t sequence_len,
        const WfGraph& graph,
        std::int32_t* score,
        const char* sequence) noexcept {
      std::uint64_t matrix_width = sequence_len + 1;
      const auto& rank_to_node = graph.rank_to_node();

      std::int32_t max_score = type_ == AlignmentType::kSW ? 0 : kNegativeInfinity;
      std::uint32_t max_i = 0;
      std::uint32_t max_j = 0;
      auto update_max_score = [&max_score, &max_i, &max_j] (
          std::int32_t* H_row,
          std::uint32_t i,
          std::uint32_t j) -> void {
        if (max_score < H_row[j]) {
          max_score = H_row[j];
          max_i = i;
          max_j = j;
        }
        return;
      };

      auto start = std::chrono::high_resolution_clock::now();
      // alignment
      for (const auto& it : rank_to_node) {
        const auto& char_profile =
            &(pimpl_->sequence_profile[it->code * matrix_width]);

        std::uint32_t i = pimpl_->node_id_to_rank[it->id] + 1;
        std::uint32_t pred_i = it->inedges.empty() ?
            0 : pimpl_->node_id_to_rank[it->inedges[0]->tail->id] + 1;

        std::int32_t* H_row = &(pimpl_->H[i * matrix_width]);
        std::int32_t* H_pred_row = &(pimpl_->H[pred_i * matrix_width]);

        // update H
        for (std::uint64_t j = 1; j < matrix_width; ++j) {
          H_row[j] = std::max(
              H_pred_row[j - 1] + char_profile[j],
              H_pred_row[j] + g_);
        }
        // check other predeccessors
        for (std::uint32_t p = 1; p < it->inedges.size(); ++p) {
          pred_i = pimpl_->node_id_to_rank[it->inedges[p]->tail->id] + 1;

          H_pred_row = &(pimpl_->H[pred_i * matrix_width]);

          for (std::uint64_t j = 1; j < matrix_width; ++j) {
            H_row[j] = std::max(
                H_pred_row[j - 1] + char_profile[j],
                std::max(
                    H_row[j],
                    H_pred_row[j] + g_));
          }
        }

        for (std::uint64_t j = 1; j < matrix_width; ++j) {
          H_row[j] = std::max(H_row[j - 1] + g_, H_row[j]);

          if (type_ == AlignmentType::kSW) {
            H_row[j] = std::max(H_row[j], 0);
            update_max_score(H_row, i, j);
          } else if (type_ == AlignmentType::kNW &&
              it->outedges.empty() && j == matrix_width - 1) {
            update_max_score(H_row, i, j);
          } 
        }
      }


      auto end = std::chrono::high_resolution_clock::now();

      // Calculate the duration
      std::chrono::duration<double, std::micro> duration = end - start;
      total_time_measure_spoa_fulltable += static_cast<double>(duration.count() );
      

      auto a2 = traceback(max_i, max_j, rank_to_node, matrix_width, pimpl_->H);


      if (mxi == 0 && mxj == 0) {
        return Alignment();
      }
      if (score) {
        *score = max_score;
      }
      
      return a2;
    }



    Alignment Debug(
        std::uint32_t sequence_len,
        const WfGraph& graph,
        std::int32_t* score,
        const char* sequence) noexcept {
      std::uint64_t matrix_width = sequence_len + 1;
      const auto& rank_to_node = graph.rank_to_node();

      std::int32_t max_score = type_ == AlignmentType::kSW ? 0 : kNegativeInfinity;
      std::uint32_t max_i = 0;
      std::uint32_t max_j = 0;
      auto update_max_score = [&max_score, &max_i, &max_j] (
          std::int32_t* H_row,
          std::uint32_t i,
          std::uint32_t j) -> void {
        if (max_score < H_row[j]) {
          max_score = H_row[j];
          max_i = i;
          max_j = j;
        }
        return;
      };

#ifdef SPOA
      auto start = std::chrono::high_resolution_clock::now();
      // alignment
      for (const auto& it : rank_to_node) {
        const auto& char_profile =
            &(pimpl_->sequence_profile[it->code * matrix_width]);

        std::uint32_t i = pimpl_->node_id_to_rank[it->id] + 1;
        std::uint32_t pred_i = it->inedges.empty() ?
            0 : pimpl_->node_id_to_rank[it->inedges[0]->tail->id] + 1;

        std::int32_t* H_row = &(pimpl_->H[i * matrix_width]);
        std::int32_t* H_pred_row = &(pimpl_->H[pred_i * matrix_width]);

        // update H
        for (std::uint64_t j = 1; j < matrix_width; ++j) {
          H_row[j] = std::max(
              H_pred_row[j - 1] + char_profile[j],
              H_pred_row[j] + g_);
        }
        // check other predeccessors
        for (std::uint32_t p = 1; p < it->inedges.size(); ++p) {
          pred_i = pimpl_->node_id_to_rank[it->inedges[p]->tail->id] + 1;

          H_pred_row = &(pimpl_->H[pred_i * matrix_width]);

          for (std::uint64_t j = 1; j < matrix_width; ++j) {
            H_row[j] = std::max(
                H_pred_row[j - 1] + char_profile[j],
                std::max(
                    H_row[j],
                    H_pred_row[j] + g_));
          }
        }

        for (std::uint64_t j = 1; j < matrix_width; ++j) {
          H_row[j] = std::max(H_row[j - 1] + g_, H_row[j]);

          if (type_ == AlignmentType::kSW) {
            H_row[j] = std::max(H_row[j], 0);
            update_max_score(H_row, i, j);
          } else if (type_ == AlignmentType::kNW &&
              it->outedges.empty() && j == matrix_width - 1) {
            update_max_score(H_row, i, j);
          } 
        }

        // print H
#ifdef debug
          std::cout << code_to_char[it->code] << "\t";
          for(int j = 0; j < matrix_width; j ++){
            std::cout << H_row[j] << "\t";
          }
          std::cout << std::endl;
#endif
      }


      auto end = std::chrono::high_resolution_clock::now();

      // Calculate the duration
      std::chrono::duration<double, std::micro> duration = end - start;


      total_time_measure_spoa_fulltable += static_cast<double>(duration.count() );
      

#endif
      std::vector<int> M(matrix_width * (rank_to_node.size() + 1), 257);

      std::string seq(sequence);
      this->seq = seq;    
      this->matrix_width = matrix_width;

      S = {};
      S.push({0, 0});
      M[0] = 0;
      d = record = 0;



      
      auto start_ = std::chrono::high_resolution_clock::now();
      while(extend(M, graph)){
          d--;
          expand(M, graph);
      } 

      auto end_ = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double, std::micro> duration_ = end_ - start_;
      total_time_measure_spoa_wavefront += static_cast<double>(duration_.count() );


      auto a2 = traceback(max_i, max_j, rank_to_node, matrix_width, pimpl_->H);
      auto al = traceback(mxi, mxj, rank_to_node, matrix_width, M.data());


#ifdef showcompare
      compare(M.data(), pimpl_->H, matrix_width, rank_to_node);
      for(auto el : a2){
          std::cout << el.first << "\t";
      }
      std::cout << std::endl;
      for(auto el : al){
          std::cout << el.first << "\t";
      }
      std::cout << std::endl << std::endl;
      for(auto el : a2){
          std::cout << el.second << "\t";
      }
      std::cout << std::endl;
      for(auto el : a2){
          std::cout << el.second << "\t";
      }
      std::cout << std::endl << std::endl;
      std::cout << max_i << ' ' << mxi << std::endl;
      std::cout << max_j << ' ' << mxj << std::endl;
#endif

      if (mxi == 0 && mxj == 0) {
        return Alignment();
      }
      if (score) {
        *score = max_score;
      }

      
      return al;
    }


    Alignment traceback(int max_i, int max_j, std::vector<WfGraph::Node*> rank_to_node, int matrix_width, int* dp){
            // backtrack
      Alignment alignment;
      std::uint32_t i = max_i;
      std::uint32_t j = max_j;
      int *tmp = pimpl_->H;
      pimpl_->H = dp;

      bool success = true;

#ifdef WFUNIT
const int G_Score = 1;
#else
const int G_Score = 2;
#endif

      auto sw_condition = [this, &i, &j, &matrix_width] () -> bool {
        return (pimpl_->H[i * matrix_width + j] == 0) ? false : true;
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
            (type_ == AlignmentType::kNW && nw_condition()) ) {
        

        auto H_ij = pimpl_->H[i * matrix_width + j];
        bool predecessor_found = false;

        if (i != 0 && j != 0) {
          const auto& it = rank_to_node[i - 1];
          std::int32_t match_cost =
              // pimpl_->sequence_profile[it->code * matrix_width + j];
              0 - G_Score * (seq[j-1] != code_to_char[it->code]);

          std::uint32_t pred_i = it->inedges.empty() ?
              0 : pimpl_->node_id_to_rank[it->inedges[0]->tail->id] + 1;

          if (H_ij == pimpl_->H[pred_i * matrix_width + (j - 1)] + match_cost) {
            prev_i = pred_i;
            prev_j = j - 1;
            predecessor_found = true;
          } else {
            for (std::uint32_t p = 1; p < it->inedges.size(); ++p) {
              std::uint32_t pred_i =
                  pimpl_->node_id_to_rank[it->inedges[p]->tail->id] + 1;

              if (H_ij == pimpl_->H[pred_i * matrix_width + (j - 1)] + match_cost) {
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

          if (H_ij == pimpl_->H[pred_i * matrix_width + j] + g_) {
            prev_i = pred_i;
            prev_j = j;
            predecessor_found = true;
          } else {
            for (std::uint32_t p = 1; p < it->inedges.size(); ++p) {
              std::uint32_t pred_i =
                  pimpl_->node_id_to_rank[it->inedges[p]->tail->id] + 1;

              if (H_ij == pimpl_->H[pred_i * matrix_width + j] + g_) {
                prev_i = pred_i;
                prev_j = j;
                predecessor_found = true;
                break;
              }
            }
          }
        }

        if (!predecessor_found && H_ij == pimpl_->H[i * matrix_width + j - 1] + g_) {  // NOLINT
          prev_i = i;
          prev_j = j - 1;
          predecessor_found = true;
        }

        alignment.emplace_back(
            i == prev_i ? -1 : rank_to_node[i - 1]->id,
            j == prev_j ? -1 : j - 1);

        if(i == prev_i && j == prev_j){
          success = false;
          // printf("debug hard q");
          break;
        }

        i = prev_i;
        j = prev_j;


      }

      std::reverse(alignment.begin(), alignment.end());

      pimpl_->H = tmp;

#ifdef showtb
        for(auto &ch : alignment){
          std::cout << ch.first << '\t';
        }
        std::cout << std::endl;

        for(auto &ch : alignment){
          std::cout << ch.second << '\t';
        }
        std::cout << std::endl << std::endl;
#endif

      if(!success){
        // printf("early return!!!!!!!!!!!!!!\n");
        Alignment al{};
        al.push_back({std::numeric_limits<int>::min(), 0});
        std::cerr << "eary return !!!!!!!!" << std::endl;
        // printf("early return!!!!!!!!!!!!!!\n");
        return al;
      }

      return alignment;
    }

    void compare(int *M, int *H, int matrix_width, std::vector<WfGraph::Node*> rank_to_node){      
      std::cout << "\t";
      for(auto c : seq){
          std::cout << "\t" << c;
      }
      std::cout << std::endl;

      for (int i = 0; i <= rank_to_node.size(); i ++){
        if(i == 0){
          std::cout << " ";
        }else{
          std::cout << code_to_char[rank_to_node[i-1]->code];
        }

        for (int j = 0; j < matrix_width; j ++){
          auto m = M[i * matrix_width + j];
          auto h = H[i * matrix_width + j];
          if(m == 257){
            std::cout << '-';
            continue;
          }
          std::cout << (m == h);
        }
        std::cout << std::endl;
      }

    }

    void print_N(int *M, std::vector<WfGraph::Node*> rank_to_node){

      std::cout << "\t";
      for(auto c : seq){
          std::cout << "\t" << c;
      }
      std::cout << std::endl;

      for (int i = 0; i <= rank_to_node.size(); i ++){
        if(i == 0){
          std::cout << "\t";
        }else{
          std::cout << code_to_char[rank_to_node[i-1]->code] << "\t";
        }

        for (int j = 0; j < matrix_width; j ++){
          auto pn = M[i * matrix_width + j];
          std::cout << ((pn >= d) ? pn : 255 ) << "\t";
        }
        std::cout << std::endl;
      }



    }
};

}  // namespace spoa

#endif  // Wfpoa_ALIGNMENT_ENGINE_HPP_
