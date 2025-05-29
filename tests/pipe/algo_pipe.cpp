#include <biovoltron/pipe/algo_pipe.hpp>
#include <catch.hpp>
#include <iostream>
#include <range/v3/action/sort.hpp>
#include <range/v3/istream_range.hpp>
#include <range/v3/view/zip.hpp>

using namespace biovoltron;

TEST_CASE("pipe") {
  std::cout << "biovoltron/pipe/algo_pipe.hpp needs test.\n";
  /*
  if (false) {
    auto fa = std::ifstream{"chr22.fa"};
    auto ref = FastaRecord<>{};
    fa >> ref;
    auto index = ref | pipe::build{.LOOKUP_LEN = 8};

    auto fq1 = std::ifstream{""};
    auto fq2 = std::ifstream{""};
    auto reads = ranges::view::zip(ranges::istream_range<FastqRecord<>>(fq1),
                                   ranges::istream_range<FastqRecord<>>(fq2));

    std::cout << ref.seq.size() << "\n";
    const auto alignments = reads
                            | pipe::align{.ref = ref,
                                          .index = std::move(index),
                                          .args = {.MAX_HIT_CNT = 512}}
                            | ranges::action::sort;

    auto sam = std::ofstream{""};
    for (const auto& alignment : alignments) sam << alignment << "\n";
    sam.close();

    auto vcf = std::ofstream{""};
    const auto variants
      = alignments
        | pipe::call{.ref = std::move(ref),
                     .args = {.MAX_READS_PER_ALIGN_BEGIN = 10}};
    for (const auto& variant : variants) vcf << variant << "\n";
  }
  */
}
