#pragma once

#include <set>
#include <map>
#include <numeric>

#include <spdlog/spdlog.h>

#include <biovoltron/file_io/sam.hpp>
#include <biovoltron/utility/haplotype/haplotype.hpp>
#include <biovoltron/utility/read/quality_utils.hpp>
#include <biovoltron/algo/assemble/graph/graph_wrapper.hpp>
#include <biovoltron/algo/align/inexact_match/smithwaterman.hpp>

namespace biovoltron {

struct HaplotypeGraph {
  struct Parameter {
    const std::size_t DEFAULT_NUM_PATHS = 128;
    const std::size_t MIN_BASE_QUALITY = 10 + QualityUtils::ASCII_OFFSET;
    const std::size_t PRUNE_FACTOR = 2;
  } para;

private:
  struct VertexProperty {
    std::string_view kmer;
  };

  struct EdgeProperty {
    int count = 0;
    bool is_ref = false;
    bool is_on_path = false;
    double score = std::numeric_limits<double>::lowest();
  };

  using Graph = GraphWrapper<VertexProperty, EdgeProperty>;
  using Vertex = Graph::Vertex;
  using Edge = Graph::Edge;
  using Path = std::vector<Vertex>;

  // TODO: Edge filter;

  Graph g;
  Vertex source{}, sink{};
  std::vector<Path> paths;
  std::set<Vertex> vertices_on_paths;

  std::string_view ref;
  std::vector<std::string_view> read_segs;

  std::size_t kmer_size;
  std::set<std::string_view> dup_kmers;
  std::map<std::string_view, Vertex> unique_kmers;

  auto create_edge(Vertex u, Vertex v, bool is_ref) {
    const auto e = g.create_edge(u, v);
    g[e].count += 1;
    g[e].is_ref = is_ref;
  }

  auto create_vertex(std::string_view kmer) {
    const auto v = g.create_vertex();
    g[v].kmer = kmer;
    if (dup_kmers.find(kmer) == dup_kmers.end()) {
      unique_kmers.emplace(kmer, v);
    }
    return v;
  }

  auto get_vertex(std::string_view kmer) {
    if (auto it = unique_kmers.find(kmer); it != unique_kmers.end()) {
      return it->second;
    }
    return create_vertex(kmer);
  }

  auto increase_counts_backwards(Vertex v, std::string_view kmer) {
    if (kmer.empty()) {
      return;
    }
    if (g.in_degree(v) == 1) {
      for (const auto e : g.in_edges(v)) {
        const auto u = g.source(e);
        if (g[u].kmer.back() == kmer.back()) {
          g[e].count++;
          increase_counts_backwards(u, kmer.substr(0, kmer.size() - 1));
        }
      }
    }
  }

  auto extend_chain(Vertex u, std::string_view kmer, bool is_ref) {
    for (const auto e : g.out_edges(u)) {
      const auto v = g.target(e);
      if (g[v].kmer.back() == kmer.back()) {
        g[e].count++;
        return v;
      }
    }
    const auto v = get_vertex(kmer);
    create_edge(u, v, is_ref);
    return v;
  }

  auto add_seq(std::string_view seq, bool is_ref) {
    auto v = get_vertex(seq.substr(0, kmer_size));
    // increase_counts_backwards(v, seq.substr(0, kmer_size-1));
    if (is_ref) {
      source = v;
    }
    for (auto i = 1; i <= seq.size() - kmer_size; i++) {
      v = extend_chain(v, seq.substr(i, kmer_size), is_ref);
    }
    if (is_ref) {
      sink = v;
    }
  }

  auto path_finder(Vertex from, Vertex to, Path& path) -> void {
    path.push_back(from);
    if (from == to) {
      paths.push_back(path);
      for (const auto v : path) { 
        vertices_on_paths.insert(v);
      }
    } else {
      for (const auto e : g.out_edges(from)) {
        if (g[e].is_ref || g[e].count >= para.PRUNE_FACTOR
            || g.out_degree(from) == 1) {
          const auto v = g.target(e);
          if (std::ranges::find(path, v) == path.end()) { 
            path_finder(v, to, path);
          }
        }
      }
    }
    path.pop_back();
  }

  auto find_all_paths() {
    auto path = Path{};
    path_finder(source, sink, path);
  }

  auto mark_edges_on_paths() {
    for (const auto& path : paths) {
      auto u = path[0];
      for (auto i = 1u; i < path.size(); i++) {
        const auto v = path[i];
        g[g.edge(u, v)].is_on_path = true;
        u = v;
      }
    }
  }

  auto compute_edges_score() {
    for (const auto v : vertices_on_paths) {
      auto edges = std::vector<Edge> {};
      for (const auto e : g.out_edges(v)) {
        if (g[e].is_on_path) {
          edges.push_back(e);
        }
      }
      auto sum = std::accumulate(edges.begin(), edges.end(), 0.0, 
        [this](double sum, const Edge& e) {
          return sum + g[e].count;
        }
      );
      std::ranges::for_each(edges, [this, &sum](const Edge& e) {
        g[e].score = std::log10(g[e].count / sum);
      });
    }
  }

  auto get_haplotypes() {
    auto haplotypes = std::vector<Haplotype>{};

    for (const auto& path : paths) {
      auto u = path[0];
      auto seq = std::string(g[u].kmer.data(), g[u].kmer.size());
      auto score = 0.0;
      for (auto i = 1; i < path.size(); i++) {
        auto v = path[i];
        seq += g[v].kmer.back();
        score += g[g.edge(u, v)].score;
        u = v;
      }
      haplotypes.push_back({.seq = std::move(seq), .score = score});
    }

    std::ranges::sort(haplotypes, std::ranges::greater {}, &Haplotype::score);

    if (haplotypes.size() > para.DEFAULT_NUM_PATHS) {
      haplotypes.resize(para.DEFAULT_NUM_PATHS);
    }
    if (haplotypes.size() > 1) {
      SPDLOG_DEBUG("Found {} candidate haplotypes.", haplotypes.size());
    } else { 
      SPDLOG_DEBUG("Found only the reference haplotype in the assembly graph.");
    }

    for (auto& h : haplotypes) {
      auto [align_begin, cigar] = SmithWaterman::align(ref, h.seq);
      // SPDLOG_DEBUG("Adding haplotype {} from graph with kmer {}", cigar, kmer_size);
      h.align_begin_wrt_ref = align_begin;
      h.cigar = std::move(cigar);
    }

    for (auto& h : haplotypes) {
      SPDLOG_DEBUG("{}", h.seq);
      // SPDLOG_DEBUG("> Cigar = {} score {}", h.cigar, h.score);
    }
    return haplotypes;
  }

public:

  HaplotypeGraph(int kmer_size) : kmer_size(kmer_size) {

  }

  static auto get_dup_kmers(std::string_view seq, int size) {
    auto all_kmers = std::set<std::string_view> {};
    auto dup_kmers = std::set<std::string_view> {};
    for (auto i = 0; i <= seq.size() - size; i++) {
      const auto kmer = seq.substr(i, size);
      if (auto [iter, success] = all_kmers.insert(kmer); !success) {
        dup_kmers.insert(kmer);
      }
    }
    return dup_kmers;
  }

  auto set_ref(std::string_view ref) {
    this->ref = ref;
  }

  auto set_read(const SamRecord<>& read) {
    auto seq = static_cast<std::string_view>(read.seq);
    const auto& qual = read.qual;

    auto start = std::string_view::npos;
    auto is_usable = [this](auto base, auto qual) {
      return base != 'N' && qual >= para.MIN_BASE_QUALITY;
    };
    for (auto i = 0; i <= seq.size(); i++) {
      if (i == seq.size() || !is_usable(seq[i], qual[i])) {
        if (start != std::string_view::npos && i - start >= kmer_size) {
          read_segs.push_back(seq.substr(start, i - start));
        }
        start = std::string_view::npos;
      } else if (start == std::string_view::npos) {
        start = i;
      }
    }
  }

  auto build() {
    for (const auto kmer : get_dup_kmers(ref, kmer_size)) {
      dup_kmers.insert(kmer);
    }
    for (const auto seg : read_segs) {
      for (const auto kmer : get_dup_kmers(seg, kmer_size)) {
        dup_kmers.insert(kmer);
      }
    }
    add_seq(ref, true);
    for (const auto seg : read_segs) { 
      add_seq(seg, false);
    }
  }

  auto has_cycles() const {
    // TODO: has_cycle
    return true;
    // auto has_cycle = false;
    // auto vis = cycle_detector {has_cycle};
    // auto fg = boost::filtered_graph<Graph, EdgeFilter>(g, filter);
    // boost::depth_first_search(fg, boost::visitor(vis));
    // return has_cycle;
  }

  auto unique_kmers_count() const {
    return unique_kmers.size();
  }

  auto find_paths() {
    find_all_paths();
    mark_edges_on_paths();
    compute_edges_score();
    return get_haplotypes();
  }
};



} // namespace biovoltron