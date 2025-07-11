#include <biovoltron/algo/assemble/assembler.hpp>
#include <biovoltron/utility/haplotype/haplotype.hpp>
#include <biovoltron/file_io/sam.hpp>
#include <catch.hpp>

using namespace biovoltron;

auto make_sam_record(const auto& seq) {
  auto record = SamRecord{};
  record.seq = seq;
  record.qual = std::string(record.seq.size(), ';');
  return record;
}

auto ngs_sequencing(const auto& seq) {

  // NGS reads not longer than 255 base 
  auto read1 = seq.substr(100, 200);
  auto read2 = seq.substr(200, 400);
  auto read3 = seq.substr(700, 900);

  auto records = std::vector<SamRecord<false>>{};
  for (auto i = 0; i < 5; i++) {
    records.push_back(make_sam_record(read1));
    records.push_back(make_sam_record(read2));
    records.push_back(make_sam_record(read3));
  }

  return records;
}

TEST_CASE("Assembler") {
  auto ref = std::string{"AATCGAAGGTCGTAAGGACACGGTTGAGCGTTCAGCGTTCATGTGAGTCCTCACCACTTATGGCTCCATAGCCTGCTATTTAAGTGGGTTACCGGTCTCCGCCAAGTAGCTGGTGTAAGAACACAGTAACTGAGCCCAGTGTGATCAGCCCTAACGAGGTACCAATCATCATAAACCATCCCTTGAGTCTCGGTTCTGCTGGTTTCGGACGTCGTGTGGATGCAGCGGATTTGACTACCGTCCTCATAGGAATGCCGGATGTAATGAAACTTCCGCTTCCAAATATACGATATCAAAGGTGTGGAGACGATGGGCGAACTTGGCAGCGGCCCCCCCCACCGGGGGTCTCCGGCGTAGCGGTACGGTCTATGCTAAGGCGTTGCTAACAATTGCAGGAGCACGGGGCTCGCAAGTAAAAGCACCGCACTGGGCATGATACCGGGGAATACGGAGTCTTCCCTTATGCCGAAAGAAGCAGCTATAACTTCCTCGGGTAAAGGGCAAGAGAAGATCGTAGGCACGTACTCCCGAACTTCAAGAGATCCCGGTTTGCTGCGCCAACCCAATGGTAGCCACATCACGCATATTAGACCGTTGCTGAAATAGTAAAGGCCGCAACCTTCAGATGTCAGCCTTTTCATGCTGTGGATTAACAAGAGTGGGGAAGCAATACGAAGTGAGTTCGTTGGGCATGCGGGAGGGCGGCAGGAAGCAAACGGGTTGCGGCCCGGCGCGGTACGTTGTGAATCGATCTCTGACGCATACCCTCCAGCAATTCCTAAAACCTCCGCATTTTTAGATGTCTGCTGTCGGTCAGGTAGTCAACAGGTTTGTTCACCGAAACGACTGGTCTTCACCCCGTCAAATCATTAAACGCGCCCGCAGTGCTTTCACGGGTCCCCGACGTCAGATCGCCCTAGACCATGATGCCCGGTACCAAAGTCTCACTGCCGTCACGGTAAGTGGTATATGCGGTTGGGCGGCTCTCTACTTCGGTTGATGAATAATGGTGCTGAAGGCGACT"};
  auto subject = ref;
  for (auto idx = 150; idx < 155; idx++) {
    ref[idx] = 'A';
    subject[idx] = 'T';
  }
  subject.insert(subject.begin() + 300, 'T');
  subject.erase(subject.begin() + 800);

  auto reads = ngs_sequencing(subject);

  auto assembler = HaplotypeAssembler{};
  auto haplotypes = assembler.assemble(reads, std::string_view(ref));

  REQUIRE(haplotypes.size() == 8);

  REQUIRE(haplotypes[0].seq == "AATCGAAGGTCGTAAGGACACGGTTGAGCGTTCAGCGTTCATGTGAGTCCTCACCACTTATGGCTCCATAGCCTGCTATTTAAGTGGGTTACCGGTCTCCGCCAAGTAGCTGGTGTAAGAACACAGTAACTGAGCCCAGTGTGATCAGCCTTTTTGAGGTACCAATCATCATAAACCATCCCTTGAGTCTCGGTTCTGCTGGTTTCGGACGTCGTGTGGATGCAGCGGATTTGACTACCGTCCTCATAGGAATGCCGGATGTAATGAAACTTCCGCTTCCAAATATACGATATCAAAGGTTGTGGAGACGATGGGCGAACTTGGCAGCGGCCCCCCCCACCGGGGGTCTCCGGCGTAGCGGTACGGTCTATGCTAAGGCGTTGCTAACAATTGCAGGAGCACGGGGCTCGCAAGTAAAAGCACCGCACTGGGCATGATACCGGGGAATACGGAGTCTTCCCTTATGCCGAAAGAAGCAGCTATAACTTCCTCGGGTAAAGGGCAAGAGAAGATCGTAGGCACGTACTCCCGAACTTCAAGAGATCCCGGTTTGCTGCGCCAACCCAATGGTAGCCACATCACGCATATTAGACCGTTGCTGAAATAGTAAAGGCCGCAACCTTCAGATGTCAGCCTTTTCATGCTGTGGATTAACAAGAGTGGGGAAGCAATACGAAGTGAGTTCGTTGGGCATGCGGGAGGGCGGCAGGAAGCAAACGGGTTGCGGCCCGGCGCGGTACGTTGTGAATCGATCTCTGACGCATACCCTCCAGCAATTCCTAAAACCTCCGCATTTTTAGTGTCTGCTGTCGGTCAGGTAGTCAACAGGTTTGTTCACCGAAACGACTGGTCTTCACCCCGTCAAATCATTAAACGCGCCCGCAGTGCTTTCACGGGTCCCCGACGTCAGATCGCCCTAGACCATGATGCCCGGTACCAAAGTCTCACTGCCGTCACGGTAAGTGGTATATGCGGTTGGGCGGCTCTCTACTTCGGTTGATGAATAATGGTGCTGAAGGCGACT");
  REQUIRE(haplotypes[0].align_begin_wrt_ref == 0);
  REQUIRE(haplotypes[0].score == Approx(-0.237543));
  REQUIRE(haplotypes[0].cigar == "150M5D5I144M1I500M1D224M");

  REQUIRE(haplotypes[1].seq == "AATCGAAGGTCGTAAGGACACGGTTGAGCGTTCAGCGTTCATGTGAGTCCTCACCACTTATGGCTCCATAGCCTGCTATTTAAGTGGGTTACCGGTCTCCGCCAAGTAGCTGGTGTAAGAACACAGTAACTGAGCCCAGTGTGATCAGCCAAAAAGAGGTACCAATCATCATAAACCATCCCTTGAGTCTCGGTTCTGCTGGTTTCGGACGTCGTGTGGATGCAGCGGATTTGACTACCGTCCTCATAGGAATGCCGGATGTAATGAAACTTCCGCTTCCAAATATACGATATCAAAGGTTGTGGAGACGATGGGCGAACTTGGCAGCGGCCCCCCCCACCGGGGGTCTCCGGCGTAGCGGTACGGTCTATGCTAAGGCGTTGCTAACAATTGCAGGAGCACGGGGCTCGCAAGTAAAAGCACCGCACTGGGCATGATACCGGGGAATACGGAGTCTTCCCTTATGCCGAAAGAAGCAGCTATAACTTCCTCGGGTAAAGGGCAAGAGAAGATCGTAGGCACGTACTCCCGAACTTCAAGAGATCCCGGTTTGCTGCGCCAACCCAATGGTAGCCACATCACGCATATTAGACCGTTGCTGAAATAGTAAAGGCCGCAACCTTCAGATGTCAGCCTTTTCATGCTGTGGATTAACAAGAGTGGGGAAGCAATACGAAGTGAGTTCGTTGGGCATGCGGGAGGGCGGCAGGAAGCAAACGGGTTGCGGCCCGGCGCGGTACGTTGTGAATCGATCTCTGACGCATACCCTCCAGCAATTCCTAAAACCTCCGCATTTTTAGTGTCTGCTGTCGGTCAGGTAGTCAACAGGTTTGTTCACCGAAACGACTGGTCTTCACCCCGTCAAATCATTAAACGCGCCCGCAGTGCTTTCACGGGTCCCCGACGTCAGATCGCCCTAGACCATGATGCCCGGTACCAAAGTCTCACTGCCGTCACGGTAAGTGGTATATGCGGTTGGGCGGCTCTCTACTTCGGTTGATGAATAATGGTGCTGAAGGCGACT");
  REQUIRE(haplotypes[1].align_begin_wrt_ref == 0);
  REQUIRE(haplotypes[1].score == Approx(-0.936514f));
  REQUIRE(haplotypes[1].cigar == "299M1I500M1D224M");

  REQUIRE(haplotypes[2].seq == "AATCGAAGGTCGTAAGGACACGGTTGAGCGTTCAGCGTTCATGTGAGTCCTCACCACTTATGGCTCCATAGCCTGCTATTTAAGTGGGTTACCGGTCTCCGCCAAGTAGCTGGTGTAAGAACACAGTAACTGAGCCCAGTGTGATCAGCCTTTTTGAGGTACCAATCATCATAAACCATCCCTTGAGTCTCGGTTCTGCTGGTTTCGGACGTCGTGTGGATGCAGCGGATTTGACTACCGTCCTCATAGGAATGCCGGATGTAATGAAACTTCCGCTTCCAAATATACGATATCAAAGGTGTGGAGACGATGGGCGAACTTGGCAGCGGCCCCCCCCACCGGGGGTCTCCGGCGTAGCGGTACGGTCTATGCTAAGGCGTTGCTAACAATTGCAGGAGCACGGGGCTCGCAAGTAAAAGCACCGCACTGGGCATGATACCGGGGAATACGGAGTCTTCCCTTATGCCGAAAGAAGCAGCTATAACTTCCTCGGGTAAAGGGCAAGAGAAGATCGTAGGCACGTACTCCCGAACTTCAAGAGATCCCGGTTTGCTGCGCCAACCCAATGGTAGCCACATCACGCATATTAGACCGTTGCTGAAATAGTAAAGGCCGCAACCTTCAGATGTCAGCCTTTTCATGCTGTGGATTAACAAGAGTGGGGAAGCAATACGAAGTGAGTTCGTTGGGCATGCGGGAGGGCGGCAGGAAGCAAACGGGTTGCGGCCCGGCGCGGTACGTTGTGAATCGATCTCTGACGCATACCCTCCAGCAATTCCTAAAACCTCCGCATTTTTAGTGTCTGCTGTCGGTCAGGTAGTCAACAGGTTTGTTCACCGAAACGACTGGTCTTCACCCCGTCAAATCATTAAACGCGCCCGCAGTGCTTTCACGGGTCCCCGACGTCAGATCGCCCTAGACCATGATGCCCGGTACCAAAGTCTCACTGCCGTCACGGTAAGTGGTATATGCGGTTGGGCGGCTCTCTACTTCGGTTGATGAATAATGGTGCTGAAGGCGACT");
  REQUIRE(haplotypes[2].align_begin_wrt_ref == 0);
  REQUIRE(haplotypes[2].score == Approx(-0.936514f));
  REQUIRE(haplotypes[2].cigar == "150M5D5I644M1D224M");

  REQUIRE(haplotypes[3].seq == "AATCGAAGGTCGTAAGGACACGGTTGAGCGTTCAGCGTTCATGTGAGTCCTCACCACTTATGGCTCCATAGCCTGCTATTTAAGTGGGTTACCGGTCTCCGCCAAGTAGCTGGTGTAAGAACACAGTAACTGAGCCCAGTGTGATCAGCCTTTTTGAGGTACCAATCATCATAAACCATCCCTTGAGTCTCGGTTCTGCTGGTTTCGGACGTCGTGTGGATGCAGCGGATTTGACTACCGTCCTCATAGGAATGCCGGATGTAATGAAACTTCCGCTTCCAAATATACGATATCAAAGGTTGTGGAGACGATGGGCGAACTTGGCAGCGGCCCCCCCCACCGGGGGTCTCCGGCGTAGCGGTACGGTCTATGCTAAGGCGTTGCTAACAATTGCAGGAGCACGGGGCTCGCAAGTAAAAGCACCGCACTGGGCATGATACCGGGGAATACGGAGTCTTCCCTTATGCCGAAAGAAGCAGCTATAACTTCCTCGGGTAAAGGGCAAGAGAAGATCGTAGGCACGTACTCCCGAACTTCAAGAGATCCCGGTTTGCTGCGCCAACCCAATGGTAGCCACATCACGCATATTAGACCGTTGCTGAAATAGTAAAGGCCGCAACCTTCAGATGTCAGCCTTTTCATGCTGTGGATTAACAAGAGTGGGGAAGCAATACGAAGTGAGTTCGTTGGGCATGCGGGAGGGCGGCAGGAAGCAAACGGGTTGCGGCCCGGCGCGGTACGTTGTGAATCGATCTCTGACGCATACCCTCCAGCAATTCCTAAAACCTCCGCATTTTTAGATGTCTGCTGTCGGTCAGGTAGTCAACAGGTTTGTTCACCGAAACGACTGGTCTTCACCCCGTCAAATCATTAAACGCGCCCGCAGTGCTTTCACGGGTCCCCGACGTCAGATCGCCCTAGACCATGATGCCCGGTACCAAAGTCTCACTGCCGTCACGGTAAGTGGTATATGCGGTTGGGCGGCTCTCTACTTCGGTTGATGAATAATGGTGCTGAAGGCGACT");
  REQUIRE(haplotypes[3].align_begin_wrt_ref == 0);
  REQUIRE(haplotypes[3].score == Approx(-0.936514f));
  REQUIRE(haplotypes[3].cigar == "150M5D5I144M1I725M");

  REQUIRE(haplotypes[4].seq == "AATCGAAGGTCGTAAGGACACGGTTGAGCGTTCAGCGTTCATGTGAGTCCTCACCACTTATGGCTCCATAGCCTGCTATTTAAGTGGGTTACCGGTCTCCGCCAAGTAGCTGGTGTAAGAACACAGTAACTGAGCCCAGTGTGATCAGCCAAAAAGAGGTACCAATCATCATAAACCATCCCTTGAGTCTCGGTTCTGCTGGTTTCGGACGTCGTGTGGATGCAGCGGATTTGACTACCGTCCTCATAGGAATGCCGGATGTAATGAAACTTCCGCTTCCAAATATACGATATCAAAGGTTGTGGAGACGATGGGCGAACTTGGCAGCGGCCCCCCCCACCGGGGGTCTCCGGCGTAGCGGTACGGTCTATGCTAAGGCGTTGCTAACAATTGCAGGAGCACGGGGCTCGCAAGTAAAAGCACCGCACTGGGCATGATACCGGGGAATACGGAGTCTTCCCTTATGCCGAAAGAAGCAGCTATAACTTCCTCGGGTAAAGGGCAAGAGAAGATCGTAGGCACGTACTCCCGAACTTCAAGAGATCCCGGTTTGCTGCGCCAACCCAATGGTAGCCACATCACGCATATTAGACCGTTGCTGAAATAGTAAAGGCCGCAACCTTCAGATGTCAGCCTTTTCATGCTGTGGATTAACAAGAGTGGGGAAGCAATACGAAGTGAGTTCGTTGGGCATGCGGGAGGGCGGCAGGAAGCAAACGGGTTGCGGCCCGGCGCGGTACGTTGTGAATCGATCTCTGACGCATACCCTCCAGCAATTCCTAAAACCTCCGCATTTTTAGATGTCTGCTGTCGGTCAGGTAGTCAACAGGTTTGTTCACCGAAACGACTGGTCTTCACCCCGTCAAATCATTAAACGCGCCCGCAGTGCTTTCACGGGTCCCCGACGTCAGATCGCCCTAGACCATGATGCCCGGTACCAAAGTCTCACTGCCGTCACGGTAAGTGGTATATGCGGTTGGGCGGCTCTCTACTTCGGTTGATGAATAATGGTGCTGAAGGCGACT");
  REQUIRE(haplotypes[4].align_begin_wrt_ref == 0);
  REQUIRE(haplotypes[4].score == Approx(-1.63548f));
  REQUIRE(haplotypes[4].cigar == "299M1I725M");

  REQUIRE(haplotypes[5].seq == "AATCGAAGGTCGTAAGGACACGGTTGAGCGTTCAGCGTTCATGTGAGTCCTCACCACTTATGGCTCCATAGCCTGCTATTTAAGTGGGTTACCGGTCTCCGCCAAGTAGCTGGTGTAAGAACACAGTAACTGAGCCCAGTGTGATCAGCCTTTTTGAGGTACCAATCATCATAAACCATCCCTTGAGTCTCGGTTCTGCTGGTTTCGGACGTCGTGTGGATGCAGCGGATTTGACTACCGTCCTCATAGGAATGCCGGATGTAATGAAACTTCCGCTTCCAAATATACGATATCAAAGGTGTGGAGACGATGGGCGAACTTGGCAGCGGCCCCCCCCACCGGGGGTCTCCGGCGTAGCGGTACGGTCTATGCTAAGGCGTTGCTAACAATTGCAGGAGCACGGGGCTCGCAAGTAAAAGCACCGCACTGGGCATGATACCGGGGAATACGGAGTCTTCCCTTATGCCGAAAGAAGCAGCTATAACTTCCTCGGGTAAAGGGCAAGAGAAGATCGTAGGCACGTACTCCCGAACTTCAAGAGATCCCGGTTTGCTGCGCCAACCCAATGGTAGCCACATCACGCATATTAGACCGTTGCTGAAATAGTAAAGGCCGCAACCTTCAGATGTCAGCCTTTTCATGCTGTGGATTAACAAGAGTGGGGAAGCAATACGAAGTGAGTTCGTTGGGCATGCGGGAGGGCGGCAGGAAGCAAACGGGTTGCGGCCCGGCGCGGTACGTTGTGAATCGATCTCTGACGCATACCCTCCAGCAATTCCTAAAACCTCCGCATTTTTAGATGTCTGCTGTCGGTCAGGTAGTCAACAGGTTTGTTCACCGAAACGACTGGTCTTCACCCCGTCAAATCATTAAACGCGCCCGCAGTGCTTTCACGGGTCCCCGACGTCAGATCGCCCTAGACCATGATGCCCGGTACCAAAGTCTCACTGCCGTCACGGTAAGTGGTATATGCGGTTGGGCGGCTCTCTACTTCGGTTGATGAATAATGGTGCTGAAGGCGACT");
  REQUIRE(haplotypes[5].align_begin_wrt_ref == 0);
  REQUIRE(haplotypes[5].score == Approx(-1.63548f));
  REQUIRE(haplotypes[5].cigar == "150M5D5I869M");

  REQUIRE(haplotypes[6].seq == "AATCGAAGGTCGTAAGGACACGGTTGAGCGTTCAGCGTTCATGTGAGTCCTCACCACTTATGGCTCCATAGCCTGCTATTTAAGTGGGTTACCGGTCTCCGCCAAGTAGCTGGTGTAAGAACACAGTAACTGAGCCCAGTGTGATCAGCCAAAAAGAGGTACCAATCATCATAAACCATCCCTTGAGTCTCGGTTCTGCTGGTTTCGGACGTCGTGTGGATGCAGCGGATTTGACTACCGTCCTCATAGGAATGCCGGATGTAATGAAACTTCCGCTTCCAAATATACGATATCAAAGGTGTGGAGACGATGGGCGAACTTGGCAGCGGCCCCCCCCACCGGGGGTCTCCGGCGTAGCGGTACGGTCTATGCTAAGGCGTTGCTAACAATTGCAGGAGCACGGGGCTCGCAAGTAAAAGCACCGCACTGGGCATGATACCGGGGAATACGGAGTCTTCCCTTATGCCGAAAGAAGCAGCTATAACTTCCTCGGGTAAAGGGCAAGAGAAGATCGTAGGCACGTACTCCCGAACTTCAAGAGATCCCGGTTTGCTGCGCCAACCCAATGGTAGCCACATCACGCATATTAGACCGTTGCTGAAATAGTAAAGGCCGCAACCTTCAGATGTCAGCCTTTTCATGCTGTGGATTAACAAGAGTGGGGAAGCAATACGAAGTGAGTTCGTTGGGCATGCGGGAGGGCGGCAGGAAGCAAACGGGTTGCGGCCCGGCGCGGTACGTTGTGAATCGATCTCTGACGCATACCCTCCAGCAATTCCTAAAACCTCCGCATTTTTAGTGTCTGCTGTCGGTCAGGTAGTCAACAGGTTTGTTCACCGAAACGACTGGTCTTCACCCCGTCAAATCATTAAACGCGCCCGCAGTGCTTTCACGGGTCCCCGACGTCAGATCGCCCTAGACCATGATGCCCGGTACCAAAGTCTCACTGCCGTCACGGTAAGTGGTATATGCGGTTGGGCGGCTCTCTACTTCGGTTGATGAATAATGGTGCTGAAGGCGACT");
  REQUIRE(haplotypes[6].align_begin_wrt_ref == 0);
  REQUIRE(haplotypes[6].score == Approx(-1.63548f));
  REQUIRE(haplotypes[6].cigar == "799M1D224M");

  REQUIRE(haplotypes[7].seq == "AATCGAAGGTCGTAAGGACACGGTTGAGCGTTCAGCGTTCATGTGAGTCCTCACCACTTATGGCTCCATAGCCTGCTATTTAAGTGGGTTACCGGTCTCCGCCAAGTAGCTGGTGTAAGAACACAGTAACTGAGCCCAGTGTGATCAGCCAAAAAGAGGTACCAATCATCATAAACCATCCCTTGAGTCTCGGTTCTGCTGGTTTCGGACGTCGTGTGGATGCAGCGGATTTGACTACCGTCCTCATAGGAATGCCGGATGTAATGAAACTTCCGCTTCCAAATATACGATATCAAAGGTGTGGAGACGATGGGCGAACTTGGCAGCGGCCCCCCCCACCGGGGGTCTCCGGCGTAGCGGTACGGTCTATGCTAAGGCGTTGCTAACAATTGCAGGAGCACGGGGCTCGCAAGTAAAAGCACCGCACTGGGCATGATACCGGGGAATACGGAGTCTTCCCTTATGCCGAAAGAAGCAGCTATAACTTCCTCGGGTAAAGGGCAAGAGAAGATCGTAGGCACGTACTCCCGAACTTCAAGAGATCCCGGTTTGCTGCGCCAACCCAATGGTAGCCACATCACGCATATTAGACCGTTGCTGAAATAGTAAAGGCCGCAACCTTCAGATGTCAGCCTTTTCATGCTGTGGATTAACAAGAGTGGGGAAGCAATACGAAGTGAGTTCGTTGGGCATGCGGGAGGGCGGCAGGAAGCAAACGGGTTGCGGCCCGGCGCGGTACGTTGTGAATCGATCTCTGACGCATACCCTCCAGCAATTCCTAAAACCTCCGCATTTTTAGATGTCTGCTGTCGGTCAGGTAGTCAACAGGTTTGTTCACCGAAACGACTGGTCTTCACCCCGTCAAATCATTAAACGCGCCCGCAGTGCTTTCACGGGTCCCCGACGTCAGATCGCCCTAGACCATGATGCCCGGTACCAAAGTCTCACTGCCGTCACGGTAAGTGGTATATGCGGTTGGGCGGCTCTCTACTTCGGTTGATGAATAATGGTGCTGAAGGCGACT");
  REQUIRE(haplotypes[7].align_begin_wrt_ref == 0);
  REQUIRE(haplotypes[7].score == Approx(-2.33445f));
  REQUIRE(haplotypes[7].cigar == "1024M");
}
