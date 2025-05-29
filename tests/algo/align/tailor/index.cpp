#include <biovoltron/algo/align/tailor/index.hpp>
#include <catch.hpp>

namespace ranges = std::ranges;
using namespace biovoltron;

TEST_CASE("Tailor index") {
  SECTION("Use std::string for sequnce") {
    auto ref = std::vector<FastaRecord<>>{
      {"chr1", "CGATCGATCGATGCATCGATAGGGGGGGG"}, // 8G at index 21
      {"chr2", "TAGGGGGGGGTTATTTTAGTGATCC"}, // 8G at index: 2
      {"chr3", "CGATTAGGGGGGGGCCGGCCGGCGCG"}, // 8G at index: 6
      {"chr4", "GGGGGGGGTCGTAGGAATAGGG"}, // discontinuous 8G
      {"chr5", "GGGGGAGCTAGTAGTACTATAC"}
    };
    auto index = Index{5};
    index.make_index(ref);
    const auto ref_size = std::accumulate(ref.begin(), ref.end(), 0,
      [](std::uint32_t init, const FastaRecord<>& record) {
        return init + record.seq.size();
    });
    REQUIRE(index.get_bwt_size() == ref_size + 1); // $ at the end.

    for (const auto& chr : ref)
      CHECK(index.get_chr_size(chr.name) == chr.seq.size());

    REQUIRE_THROWS_WITH(index.get_chr_size("gg"), "Chromosome is not in the index.");

    const auto read = Codec::to_istring("GGGGGGGG");
    const auto [begin, end, offset] = index.get_range(read, 0);
    auto hits = index.get_intervals(begin, end, read.size());
    REQUIRE(hits.size() == 4);
    CHECK(ranges::find(hits, Interval{"chr1", 21, 29}) != hits.end());
    CHECK(ranges::find(hits, Interval{"chr2", 2, 10}) != hits.end());
    CHECK(ranges::find(hits, Interval{"chr3", 6, 14}) != hits.end());
    CHECK(ranges::find(hits, Interval{"chr4", 0, 8}) != hits.end());
  }

  SECTION("Use istring for sequnce") {
    auto ref = std::vector<FastaRecord<true>>{
      {"chr1", Codec::to_istring("CGATCGATCGATGCATCGATAGGGGGGGG")}, // 8G at index 21
      {"chr2", Codec::to_istring("TAGGGGGGGGTTATTTTAGTGATCC")}, // 8G at index: 2
      {"chr3", Codec::to_istring("CGATTAGGGGGGGGCCGGCCGGCGCG")}, // 8G at index: 6
      {"chr4", Codec::to_istring("GGGGGGGGTCGTAGGAATAGGG")}, // discontinuous 8G
      {"chr5", Codec::to_istring("GGGGGAGCTAGTAGTACTATAC")}
    };
    auto index = Index{5};
    index.make_index(ref);
    const auto ref_size = std::accumulate(ref.begin(), ref.end(), 0,
      [](std::uint32_t init, const FastaRecord<true>& record) {
        return init + record.seq.size();
    });
    REQUIRE(index.get_bwt_size() == ref_size + 1); // $ at the end.

    for (const auto& chr : ref)
      CHECK(index.get_chr_size(chr.name) == chr.seq.size());

    REQUIRE_THROWS_WITH(index.get_chr_size("gg"), "Chromosome is not in the index.");

    const auto read = Codec::to_istring("GGGGGGGG");
    const auto [begin, end, offset] = index.get_range(read, 0);
    auto hits = index.get_intervals(begin, end, read.size());
    REQUIRE(hits.size() == 4);
    CHECK(ranges::find(hits, Interval{"chr1", 21, 29}) != hits.end());
    CHECK(ranges::find(hits, Interval{"chr2", 2, 10}) != hits.end());
    CHECK(ranges::find(hits, Interval{"chr3", 6, 14}) != hits.end());
    CHECK(ranges::find(hits, Interval{"chr4", 0, 8}) != hits.end());
  }

  SECTION("Save/load index") {
    auto ref = std::vector<FastaRecord<true>>{
      {"chr1", Codec::to_istring("CGATCGATCGATGCATCGATAGGGGGGGG")}, // 8G at index 21
      {"chr2", Codec::to_istring("TAGGGGGGGGTTATTTTAGTGATCC")}, // 8G at index: 2
      {"chr3", Codec::to_istring("CGATTAGGGGGGGGCCGGCCGGCGCG")}, // 8G at index: 6
      {"chr4", Codec::to_istring("GGGGGGGGTCGTAGGAATAGGG")}, // discontinuous 8G
      {"chr5", Codec::to_istring("GGGGGAGCTAGTAGTACTATAC")}
    };

    auto index = Index{5};
    index.make_index(ref);
    {
      auto ofs = std::ofstream{"tailor.idx"};
      index.save(ofs);
    }

    auto index2 = Index{};
    {
      auto ifs = std::ifstream{"tailor.idx"};
      auto ofs = std::ofstream{"tailor2.idx"};
      index2.load(ifs);
      index2.save(ofs);
    }

    auto ifs1 = std::ifstream{"tailor.idx"};
    auto ifs2 = std::ifstream{"tailor2.idx"};
    auto line1 = std::string{};
    auto line2 = std::string{};
    while (std::getline(ifs1, line1) && std::getline(ifs2, line2))
      REQUIRE(line1 == line2);

    for (const auto& chr : ref)
      CHECK(index.get_chr_size(chr.name) == chr.seq.size());

    REQUIRE_THROWS_WITH(index.get_chr_size("gg"), "Chromosome is not in the index.");
  }
}
