#include <catch.hpp>
#include <filesystem>
#include <fstream>
#include <ranges>
#include <range/v3/all.hpp>
#include <iostream>
#include <biovoltron/utility/archive/gzstream.hpp>
#include <biovoltron/file_io/fastq.hpp>

using namespace biovoltron;

const auto data_path = std::filesystem::path{DATA_PATH};

TEST_CASE("gzstream") {
  SECTION("igzstream") {
    auto igs = igzstream(data_path / "test1.fastq.gz");
    auto is = std::ifstream(data_path / "test1.fastq");

    auto fqs = ranges::view::zip(
                 ranges::istream_range<FastqRecord<>>(is),
                 ranges::istream_range<FastqRecord<>>(igs)
               );

    for (auto [fq, fqgz] : fqs) {
      REQUIRE(fq.name == fqgz.name);
      REQUIRE(fq.seq == fqgz.seq);
      REQUIRE(fq.qual == fqgz.qual);
    }
  }

  SECTION("ogzstream") {
    auto is = std::ifstream(data_path / "test1.fastq");
    auto ogs = ogzstream(data_path / "test1.fastq.gzstream");
    bool f = 1;
    for (auto &fq : ranges::istream_range<FastqRecord<>>(is)) {
      if (f)
        f = 0;
      else
        ogs << '\n';
      ogs << fq;
    }

    auto igs1 = igzstream(data_path / "test1.fastq.gz");
    auto igs2 = igzstream(data_path / "test1.fastq.gzstream");

    auto fqs = ranges::view::zip(
                 ranges::istream_range<FastqRecord<>>(igs1),
                 ranges::istream_range<FastqRecord<>>(igs2)
               );

    for (auto [fq1, fq2] : fqs) {
      REQUIRE(fq1.name == fq2.name);
      REQUIRE(fq1.seq  == fq2.seq);
      REQUIRE(fq1.qual == fq2.qual);
    }
  }
}

