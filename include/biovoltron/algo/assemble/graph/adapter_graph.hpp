#pragma once

#include <biovoltron/algo/assemble/graph/graph_wrapper.hpp>
#include <biovoltron/file_io/all.hpp>
#include <spdlog/spdlog.h>


namespace biovoltron {

template<std::ranges::view SeqView>
struct AdapterGraph {
private:
  struct VertexProperty {
    SeqView kmer;
  };

  struct EdgeProperty {
    std::size_t count = 0;
  };

  using Graph = GraphWrapper<VertexProperty, EdgeProperty>;
  using Vertex = Graph::Vertex;
  using Edge = Graph::Edge;
  using Path = std::vector<Vertex>;
  using StrType = std::conditional_t<
    std::same_as<SeqView, istring_view>,
    istring,
    std::string
  >;

  Graph graph;
  std::map<SeqView, Vertex> unique_kmers;
  std::set<SeqView> dup_kmers;
  std::size_t kmer_size;
  std::size_t minimum_occurance;

  auto build_dup_kmers(const SeqView& seq) {
    std::set<SeqView> kmers;
    for (auto i = 0u; i + kmer_size <= seq.size(); i++) {
      auto subseq = seq.substr(i, kmer_size);
      if (auto it = kmers.find(subseq); it != kmers.end()) {
        dup_kmers.insert(subseq);
      }
      kmers.insert(subseq);
    }
  }

  auto create_vertex(const SeqView& kmer) {
    auto v = graph.create_vertex();
    graph[v].kmer = kmer;
    if (dup_kmers.find(kmer) == dup_kmers.end()) {
      unique_kmers[kmer] = v;
    }
    return v;
  }

  auto create_edge(const Vertex& u, const Vertex& v) {
    auto e = graph.create_edge(u, v);
    graph[e].count = 1;
    return e;
  };

  auto get_vertex(const SeqView& kmer) {
    /* if the vertex have created before */
    if (auto it = unique_kmers.find(kmer); it != unique_kmers.end()) {
      return it->second;
    }
    return create_vertex(kmer);
  };

  /**
   * [New]
   */
  auto concat_vertices(const Path& path) {
    auto u = path[0];
    auto seq = StrType { graph[u].kmer };
    for (auto i = 1u; i < path.size(); i++) {
      auto v = path[i];
      seq += graph[v].kmer.back();
    }
    return seq;
  }


  /**
   * @brief extend the current kmer to next kmer, for example, if current kmer =
   * "ACTGC", then next kmer must be "CTGC" + [ACTG], if next kmer already 
   * existed and connected, then increase the edge between current kmer and next
   * kmer by 1, otherwise, create a new edge and set its count to 1
   * @param u 
   * @param kmer 
   * @return Vertex 
   */
  auto extend_chain(const Vertex& u, const SeqView kmer) {
    for (auto e : graph.out_edges(u, false)) {
      auto v = graph.target(e);
      if (graph[v].kmer.back() == kmer.back()) {
        graph[e].count += 1;
        return v;
      }
    }
    auto v = get_vertex(kmer);
    create_edge(u, v);
    return v;
  }

  void add_seq(const SeqView seq) {
    auto v = get_vertex(seq.substr(0, kmer_size));
    for (auto i = 1u; i + kmer_size <= seq.size(); i++) {
      v = extend_chain(v, seq.substr(i, kmer_size));
    }
  }

public:
  AdapterGraph(const std::size_t kmer_size,
               const std::size_t minimum_occurance) {
    this->kmer_size = kmer_size;
    this->minimum_occurance = minimum_occurance;

    graph.set_edge_filter([this](const Edge& e) {
      return graph[e].count >= this->minimum_occurance;
    });
  }

  auto build(const std::vector<SeqView>& seqs) {
    for (const auto& seq : seqs) {
      build_dup_kmers(seq);
    }
    for (const auto& seq : seqs) {
      if (seq.size() >= kmer_size) {
        add_seq(seq);
      }
    }
  }

  auto get_adapters() -> std::vector<StrType> {
    auto all_paths = std::vector<Path>{};
    for (auto u : graph.get_sources()) {
      for (auto v : graph.get_sinks()) {
        auto paths = graph.find_paths(u, v);
        std::ranges::move(paths, std::back_inserter(all_paths));
      }
    }
    if (all_paths.empty()) {
      return {};
    }
    // get the highest frequency adapter
    std::ranges::sort(all_paths, [this](const Path &lhs, const Path& rhs) {
      // ! weird sort policy, need check
      for (const auto lhs_e : graph.out_edges(lhs[0])) {
        for (const auto rhs_e : graph.out_edges(rhs[0])) {
          // ! should consider the situation that count is same
          return graph[lhs_e].count > graph[rhs_e].count;
        }
      }
      return lhs.size() > rhs.size();
    });
    auto adapters = std::vector<StrType>{};
    for (const auto& path : all_paths) {
      adapters.emplace_back(concat_vertices(path));
    }
    return adapters;
  }
};

} // namespace biovoltron