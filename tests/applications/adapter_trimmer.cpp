#include <biovoltron/applications/adapter_trimmer/adapter_trimmer.hpp>
#include <biovoltron/file_io/all.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/ostream_sink.h>
#include <catch.hpp>
#include <iostream>
#include <fstream>
#include <ranges>
#include <algorithm>

using namespace biovoltron;

const auto data_path = std::filesystem::path{ DATA_PATH } / "adapter_trimmer";
using Fastq = FastqRecord<>;

std::ostringstream oss;
auto ostream_sink = std::make_shared<spdlog::sinks::ostream_sink_mt>(oss);
auto logger = std::make_shared<spdlog::logger>("my_logger", ostream_sink);

auto read_fastq = [](std::filesystem::path p) {
  auto fin = std::ifstream(p);
  REQUIRE(fin.is_open());
  std::vector<Fastq> v;
  Fastq fq;
  while (fin >> fq) {
    v.emplace_back(std::move(fq));
  }
  return v;
};

TEST_CASE("PairedEndAdapterTrimmer::trim - Trims adapters from paired-end reads", "[PairedEndAdapterTrimmer]") {

  spdlog::set_default_logger(logger);
  spdlog::set_level(spdlog::level::debug);

  SECTION("Normal Mode: Has adapter") {
    auto fwd_reads_path = data_path / "has_adapter_1.fq";
    auto rev_reads_path = data_path / "has_adapter_2.fq";
    auto fwd_ans_path = data_path / "trimmed_pe_1.fastq";
    auto rev_ans_path = data_path / "trimmed_pe_2.fastq";

    auto fwd_reads = read_fastq(fwd_reads_path);
    auto rev_reads = read_fastq(rev_reads_path);
    auto fwd_ans = read_fastq(fwd_ans_path);
    auto rev_ans = read_fastq(rev_ans_path);

   PairedEndAdapterTrimmer<Fastq> trimmer;
    auto&& [fwd_trim_reads, rev_trim_reads] = trimmer.trim(fwd_reads, rev_reads);

    REQUIRE(fwd_trim_reads.size() == fwd_reads.size());
    REQUIRE(fwd_trim_reads.size() == rev_trim_reads.size());
    
    REQUIRE(std::ranges::equal(fwd_ans, fwd_trim_reads, [&](auto &a, auto &b) {
      return a.seq == b.seq && a.name == b.name;
    }));
    REQUIRE(std::ranges::equal(rev_ans, rev_trim_reads, [&](auto& a, auto& b) {
      return a.seq == b.seq && a.name == b.name;
    }));
  }
  
  SECTION("Normal Mode: Has no adapter") {
    auto fwd_reads_path = data_path / "no_adapter_1.fq";
    auto rev_reads_path = data_path / "no_adapter_2.fq";

    auto fwd_reads = read_fastq(fwd_reads_path);
    auto rev_reads = read_fastq(rev_reads_path);

    PairedEndAdapterTrimmer<Fastq> trimmer;
    auto&& [fwd_trim_reads, rev_trim_reads] = trimmer.trim(fwd_reads, rev_reads);

    REQUIRE(fwd_trim_reads.size() == fwd_reads.size());
    REQUIRE(fwd_trim_reads.size() == rev_trim_reads.size());
    REQUIRE(std::ranges::equal(fwd_reads, fwd_trim_reads, [&](auto& a, auto& b) {
      return a.seq == b.seq && a.name == b.name;
    }));
    REQUIRE(std::ranges::equal(rev_reads, rev_trim_reads, [&](auto& a, auto& b) {
      return a.seq == b.seq && a.name == b.name;
    }));
  }

  SECTION("Asio mode: Has adapter") {
    auto fwd_reads_path = data_path / "has_adapter_1.fq";
    auto rev_reads_path = data_path / "has_adapter_2.fq";
    auto fwd_output_path = data_path / "has_adapter_output.fq";
    auto rev_output_path = data_path / "has_adapter_putput.fq";

    biovoltron::PairedEndAdapterTrimmer<Fastq> trimmer;
    trimmer.trim(fwd_reads_path, rev_reads_path, fwd_output_path, rev_output_path);

    auto fwd_ans_path = data_path / "trimmed_pe_1.fastq";
    auto rev_ans_path = data_path / "trimmed_pe_2.fastq";

    auto fwd_ref_reads = read_fastq(fwd_ans_path);
    auto rev_ref_reads = read_fastq(rev_ans_path);
    auto fwd_trim_reads = read_fastq(fwd_output_path);
    auto rev_trim_reads = read_fastq(rev_output_path);

    REQUIRE(fwd_ref_reads.size() == fwd_trim_reads.size());
    REQUIRE(fwd_trim_reads.size() == rev_trim_reads.size());

    REQUIRE(std::ranges::equal(fwd_ref_reads, fwd_trim_reads, [&](auto& a, auto& b) {
      return a.seq == b.seq && a.name == b.name;
    }));
    REQUIRE(std::ranges::equal(rev_ref_reads, rev_trim_reads, [&](auto& a, auto& b) {
      return a.seq == b.seq && a.name == b.name;
    }));
    std::filesystem::remove(fwd_output_path);
    std::filesystem::remove(rev_output_path);
  }

  SECTION("Asio mode: has no adapter") {
    auto fwd_reads_path = data_path / "no_adapter_1.fq";
    auto rev_reads_path = data_path / "no_adapter_2.fq";
    auto fwd_output_path = data_path / "has_adapter_output.fq";
    auto rev_output_path = data_path / "has_adapter_putput.fq";

    biovoltron::PairedEndAdapterTrimmer<Fastq> trimmer;
    trimmer.trim(fwd_reads_path, rev_reads_path, fwd_output_path, rev_output_path);

    auto fwd_ref_reads = read_fastq(fwd_reads_path);
    auto rev_ref_reads = read_fastq(rev_reads_path);
    auto fwd_trim_reads = read_fastq(fwd_output_path);
    auto rev_trim_reads = read_fastq(rev_output_path);

    REQUIRE(std::ranges::equal(fwd_ref_reads, fwd_trim_reads, [&](auto& a, auto& b) {
      return a.seq == b.seq && a.name == b.name;
    }));
    REQUIRE(std::ranges::equal(rev_ref_reads, rev_trim_reads, [&](auto& a, auto& b) {
      return a.seq == b.seq && a.name == b.name;
    }));
    std::filesystem::remove(fwd_output_path);
    std::filesystem::remove(rev_output_path);
  }
}


TEST_CASE("SingleEndAdapterTrimmer::trim - Trims adapters from single-end reads", "[SingleEndAdapterTrimmer]") {
  spdlog::set_default_logger(logger);
  spdlog::set_level(spdlog::level::debug);

  SECTION("has adapter") {
    SingleEndAdapterTrimmer<FastqRecord<>> trimmer;
    
  }
}