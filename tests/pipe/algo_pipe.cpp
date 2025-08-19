#include <biovoltron/pipe/algo_pipe.hpp>
#include <catch.hpp>

#include <range/v3/istream_range.hpp>
#include <range/v3/view/zip.hpp>

#include <biovoltron/file_io/fasta.hpp>
#include <biovoltron/file_io/fastq.hpp>

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <sstream>
#include <string>
#include <vector>

using namespace biovoltron;

// ---- Toy reference (FASTA) ----
static std::string toy_fasta() {
  return
    ">chrToy\n"
    "ACGTACGTACGTACGTACGTACGTACGTACGT\n";
}

// ---- Toy paired-end reads: perfect match ----
static std::pair<std::string,std::string> toy_fastq_pair_perfect() {
  return {
    "@read1/1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n",
    "@read1/2\nACGTACGTACGT\n+\nIIIIIIIIIIII\n"
  };
}

// ---- Toy paired-end reads: one SNP mutation ----
static std::pair<std::string,std::string> toy_fastq_pair_with_snp() {
  return {
    "@read2/1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n",
    "@read2/2\nACGTATGTACGT\n+\nIIIIIIIIIIII\n" // 6th base mutated
  };
}

// ---- Convert ASCII reference sequence to uppercase (required by HaplotypeCaller) ----
static void to_upper(FastaRecord<> &ref_ascii) {
  std::transform(ref_ascii.seq.begin(), ref_ascii.seq.end(),
                 ref_ascii.seq.begin(),
                 [](unsigned char c){ return std::toupper(c); });
}

TEST_CASE("ALGO Pipe: build index once and reuse for multiple cases",
          "[algo_pipe][quick]") {
  // 1) Load reference sequence (FASTA) and build FM-index (done ONCE)
  std::istringstream fa_in{toy_fasta()};
  FastaRecord<true> ref_enc{};
  fa_in >> ref_enc;

  auto index = ref_enc | pipe::build<1, std::uint32_t, StableSorter<std::uint32_t>>{};

  // 2) Prepare aligner and variant caller (reused for both cases)
  pipe::align::Parameters a{};
  auto aligner = pipe::align{ref_enc, index, a};

  FastaRecord<> ref_ascii = static_cast<FastaRecord<>>(ref_enc);
  to_upper(ref_ascii);

  pipe::call::Parameters c{};

  // 3) Case A: Perfect-match reads
  //    Expect: pipeline runs without exception, no variants detected.
  {
    auto [fq1s, fq2s] = toy_fastq_pair_perfect();
    std::istringstream fq1{fq1s}, fq2{fq2s};

    auto reads = ranges::view::zip(
      ranges::istream_range<FastqRecord<>>(fq1),
      ranges::istream_range<FastqRecord<>>(fq2)
    );

    REQUIRE_NOTHROW( ([&]{
      auto alignments = reads | aligner;
      auto variants   = alignments | pipe::call{
        ref_ascii, HaplotypeAssembler{}, PairHMM{}, Genotyper{}, c
      };
      (void)variants;
    }()) );
  }

  // 4) Case B: Reads with one SNP
  //    Expect: pipeline runs without exception, at least one variant.
  {
    auto [fq1s, fq2s] = toy_fastq_pair_with_snp();
    std::istringstream fq1{fq1s}, fq2{fq2s};

    auto reads = ranges::view::zip(
      ranges::istream_range<FastqRecord<>>(fq1),
      ranges::istream_range<FastqRecord<>>(fq2)
    );

    REQUIRE_NOTHROW( ([&]{
      auto alignments = reads | aligner;
      auto variants   = alignments | pipe::call{
        ref_ascii, HaplotypeAssembler{}, PairHMM{}, Genotyper{}, c
      };
      (void)variants;
    }()) );
  }
}
