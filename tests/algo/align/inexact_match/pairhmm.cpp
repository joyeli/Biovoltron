#include <biovoltron/algo/align/inexact_match/pairhmm.hpp>
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

TEST_CASE("Pairhmm") {
  auto ref = std::string{"AATCGAAGGTCGTAAGGACACGGTTGAGCGTTCAGCGTTCATGTGAGTCCTCACCACTTATGGCTCCATAGCCTGCTATTTAAGTGGGTTACCGGTCTCCGCCAAGTAGCTGGTGTAAGAACACAGTAACTGAGCCCAGTGTGATCAGCCCTAACGAGGTACCAATCATCATAAACCATCCCTTGAGTCTCGGTTCTGCTGGTTTCGGACGTCGTGTGGATGCAGCGGATTTGACTACCGTCCTCATAGGAATGCCGGATGTAATGAAACTTCCGCTTCCAAATATACGATATCAAAGGTGTGGAGACGATGGGCGAACTTGGCAGCGGCCCCCCCCACCGGGGGTCTCCGGCGTAGCGGTACGGTCTATGCTAAGGCGTTGCTAACAATTGCAGGAGCACGGGGCTCGCAAGTAAAAGCACCGCACTGGGCATGATACCGGGGAATACGGAGTCTTCCCTTATGCCGAAAGAAGCAGCTATAACTTCCTCGGGTAAAGGGCAAGAGAAGATCGTAGGCACGTACTCCCGAACTTCAAGAGATCCCGGTTTGCTGCGCCAACCCAATGGTAGCCACATCACGCATATTAGACCGTTGCTGAAATAGTAAAGGCCGCAACCTTCAGATGTCAGCCTTTTCATGCTGTGGATTAACAAGAGTGGGGAAGCAATACGAAGTGAGTTCGTTGGGCATGCGGGAGGGCGGCAGGAAGCAAACGGGTTGCGGCCCGGCGCGGTACGTTGTGAATCGATCTCTGACGCATACCCTCCAGCAATTCCTAAAACCTCCGCATTTTTAGATGTCTGCTGTCGGTCAGGTAGTCAACAGGTTTGTTCACCGAAACGACTGGTCTTCACCCCGTCAAATCATTAAACGCGCCCGCAGTGCTTTCACGGGTCCCCGACGTCAGATCGCCCTAGACCATGATGCCCGGTACCAAAGTCTCACTGCCGTCACGGTAAGTGGTATATGCGGTTGGGCGGCTCTCTACTTCGGTTGATGAATAATGGTGCTGAAGGCGACT"};
  auto subject = ref;
  for (auto idx = 150; idx < 155; idx++) {
    ref[idx] = 'A';
    subject[idx] = 'T';
  }
  subject.insert(subject.begin() + 300, 'T');
  subject.erase(subject.begin() + 800);

  auto reads = ngs_sequencing(subject);

  std::vector<Haplotype> haplotypes(8);
  haplotypes[0].seq = "AATCGAAGGTCGTAAGGACACGGTTGAGCGTTCAGCGTTCATGTGAGTCCTCACCACTTATGGCTCCATAGCCTGCTATTTAAGTGGGTTACCGGTCTCCGCCAAGTAGCTGGTGTAAGAACACAGTAACTGAGCCCAGTGTGATCAGCCTTTTTGAGGTACCAATCATCATAAACCATCCCTTGAGTCTCGGTTCTGCTGGTTTCGGACGTCGTGTGGATGCAGCGGATTTGACTACCGTCCTCATAGGAATGCCGGATGTAATGAAACTTCCGCTTCCAAATATACGATATCAAAGGTTGTGGAGACGATGGGCGAACTTGGCAGCGGCCCCCCCCACCGGGGGTCTCCGGCGTAGCGGTACGGTCTATGCTAAGGCGTTGCTAACAATTGCAGGAGCACGGGGCTCGCAAGTAAAAGCACCGCACTGGGCATGATACCGGGGAATACGGAGTCTTCCCTTATGCCGAAAGAAGCAGCTATAACTTCCTCGGGTAAAGGGCAAGAGAAGATCGTAGGCACGTACTCCCGAACTTCAAGAGATCCCGGTTTGCTGCGCCAACCCAATGGTAGCCACATCACGCATATTAGACCGTTGCTGAAATAGTAAAGGCCGCAACCTTCAGATGTCAGCCTTTTCATGCTGTGGATTAACAAGAGTGGGGAAGCAATACGAAGTGAGTTCGTTGGGCATGCGGGAGGGCGGCAGGAAGCAAACGGGTTGCGGCCCGGCGCGGTACGTTGTGAATCGATCTCTGACGCATACCCTCCAGCAATTCCTAAAACCTCCGCATTTTTAGTGTCTGCTGTCGGTCAGGTAGTCAACAGGTTTGTTCACCGAAACGACTGGTCTTCACCCCGTCAAATCATTAAACGCGCCCGCAGTGCTTTCACGGGTCCCCGACGTCAGATCGCCCTAGACCATGATGCCCGGTACCAAAGTCTCACTGCCGTCACGGTAAGTGGTATATGCGGTTGGGCGGCTCTCTACTTCGGTTGATGAATAATGGTGCTGAAGGCGACT";
  haplotypes[0].align_begin_wrt_ref = 0;
  haplotypes[0].score = -0.237544;
  haplotypes[0].cigar = "150M5D5I144M1I500M1D224M";

  haplotypes[1].seq = "AATCGAAGGTCGTAAGGACACGGTTGAGCGTTCAGCGTTCATGTGAGTCCTCACCACTTATGGCTCCATAGCCTGCTATTTAAGTGGGTTACCGGTCTCCGCCAAGTAGCTGGTGTAAGAACACAGTAACTGAGCCCAGTGTGATCAGCCAAAAAGAGGTACCAATCATCATAAACCATCCCTTGAGTCTCGGTTCTGCTGGTTTCGGACGTCGTGTGGATGCAGCGGATTTGACTACCGTCCTCATAGGAATGCCGGATGTAATGAAACTTCCGCTTCCAAATATACGATATCAAAGGTTGTGGAGACGATGGGCGAACTTGGCAGCGGCCCCCCCCACCGGGGGTCTCCGGCGTAGCGGTACGGTCTATGCTAAGGCGTTGCTAACAATTGCAGGAGCACGGGGCTCGCAAGTAAAAGCACCGCACTGGGCATGATACCGGGGAATACGGAGTCTTCCCTTATGCCGAAAGAAGCAGCTATAACTTCCTCGGGTAAAGGGCAAGAGAAGATCGTAGGCACGTACTCCCGAACTTCAAGAGATCCCGGTTTGCTGCGCCAACCCAATGGTAGCCACATCACGCATATTAGACCGTTGCTGAAATAGTAAAGGCCGCAACCTTCAGATGTCAGCCTTTTCATGCTGTGGATTAACAAGAGTGGGGAAGCAATACGAAGTGAGTTCGTTGGGCATGCGGGAGGGCGGCAGGAAGCAAACGGGTTGCGGCCCGGCGCGGTACGTTGTGAATCGATCTCTGACGCATACCCTCCAGCAATTCCTAAAACCTCCGCATTTTTAGTGTCTGCTGTCGGTCAGGTAGTCAACAGGTTTGTTCACCGAAACGACTGGTCTTCACCCCGTCAAATCATTAAACGCGCCCGCAGTGCTTTCACGGGTCCCCGACGTCAGATCGCCCTAGACCATGATGCCCGGTACCAAAGTCTCACTGCCGTCACGGTAAGTGGTATATGCGGTTGGGCGGCTCTCTACTTCGGTTGATGAATAATGGTGCTGAAGGCGACT";
  haplotypes[1].align_begin_wrt_ref = 0;
  haplotypes[1].score = -0.936514;
  haplotypes[1].cigar = "299M1I500M1D224M";

  haplotypes[2].seq = "AATCGAAGGTCGTAAGGACACGGTTGAGCGTTCAGCGTTCATGTGAGTCCTCACCACTTATGGCTCCATAGCCTGCTATTTAAGTGGGTTACCGGTCTCCGCCAAGTAGCTGGTGTAAGAACACAGTAACTGAGCCCAGTGTGATCAGCCTTTTTGAGGTACCAATCATCATAAACCATCCCTTGAGTCTCGGTTCTGCTGGTTTCGGACGTCGTGTGGATGCAGCGGATTTGACTACCGTCCTCATAGGAATGCCGGATGTAATGAAACTTCCGCTTCCAAATATACGATATCAAAGGTGTGGAGACGATGGGCGAACTTGGCAGCGGCCCCCCCCACCGGGGGTCTCCGGCGTAGCGGTACGGTCTATGCTAAGGCGTTGCTAACAATTGCAGGAGCACGGGGCTCGCAAGTAAAAGCACCGCACTGGGCATGATACCGGGGAATACGGAGTCTTCCCTTATGCCGAAAGAAGCAGCTATAACTTCCTCGGGTAAAGGGCAAGAGAAGATCGTAGGCACGTACTCCCGAACTTCAAGAGATCCCGGTTTGCTGCGCCAACCCAATGGTAGCCACATCACGCATATTAGACCGTTGCTGAAATAGTAAAGGCCGCAACCTTCAGATGTCAGCCTTTTCATGCTGTGGATTAACAAGAGTGGGGAAGCAATACGAAGTGAGTTCGTTGGGCATGCGGGAGGGCGGCAGGAAGCAAACGGGTTGCGGCCCGGCGCGGTACGTTGTGAATCGATCTCTGACGCATACCCTCCAGCAATTCCTAAAACCTCCGCATTTTTAGTGTCTGCTGTCGGTCAGGTAGTCAACAGGTTTGTTCACCGAAACGACTGGTCTTCACCCCGTCAAATCATTAAACGCGCCCGCAGTGCTTTCACGGGTCCCCGACGTCAGATCGCCCTAGACCATGATGCCCGGTACCAAAGTCTCACTGCCGTCACGGTAAGTGGTATATGCGGTTGGGCGGCTCTCTACTTCGGTTGATGAATAATGGTGCTGAAGGCGACT";
  haplotypes[2].align_begin_wrt_ref = 0;
  haplotypes[2].score = -0.936514;
  haplotypes[2].cigar = "150M5D5I644M1D224M";

  haplotypes[3].seq = "AATCGAAGGTCGTAAGGACACGGTTGAGCGTTCAGCGTTCATGTGAGTCCTCACCACTTATGGCTCCATAGCCTGCTATTTAAGTGGGTTACCGGTCTCCGCCAAGTAGCTGGTGTAAGAACACAGTAACTGAGCCCAGTGTGATCAGCCTTTTTGAGGTACCAATCATCATAAACCATCCCTTGAGTCTCGGTTCTGCTGGTTTCGGACGTCGTGTGGATGCAGCGGATTTGACTACCGTCCTCATAGGAATGCCGGATGTAATGAAACTTCCGCTTCCAAATATACGATATCAAAGGTTGTGGAGACGATGGGCGAACTTGGCAGCGGCCCCCCCCACCGGGGGTCTCCGGCGTAGCGGTACGGTCTATGCTAAGGCGTTGCTAACAATTGCAGGAGCACGGGGCTCGCAAGTAAAAGCACCGCACTGGGCATGATACCGGGGAATACGGAGTCTTCCCTTATGCCGAAAGAAGCAGCTATAACTTCCTCGGGTAAAGGGCAAGAGAAGATCGTAGGCACGTACTCCCGAACTTCAAGAGATCCCGGTTTGCTGCGCCAACCCAATGGTAGCCACATCACGCATATTAGACCGTTGCTGAAATAGTAAAGGCCGCAACCTTCAGATGTCAGCCTTTTCATGCTGTGGATTAACAAGAGTGGGGAAGCAATACGAAGTGAGTTCGTTGGGCATGCGGGAGGGCGGCAGGAAGCAAACGGGTTGCGGCCCGGCGCGGTACGTTGTGAATCGATCTCTGACGCATACCCTCCAGCAATTCCTAAAACCTCCGCATTTTTAGATGTCTGCTGTCGGTCAGGTAGTCAACAGGTTTGTTCACCGAAACGACTGGTCTTCACCCCGTCAAATCATTAAACGCGCCCGCAGTGCTTTCACGGGTCCCCGACGTCAGATCGCCCTAGACCATGATGCCCGGTACCAAAGTCTCACTGCCGTCACGGTAAGTGGTATATGCGGTTGGGCGGCTCTCTACTTCGGTTGATGAATAATGGTGCTGAAGGCGACT";
  haplotypes[3].align_begin_wrt_ref = 0;
  haplotypes[3].score = -0.936514;
  haplotypes[3].cigar = "150M5D5I144M1I725M";

  haplotypes[4].seq = "AATCGAAGGTCGTAAGGACACGGTTGAGCGTTCAGCGTTCATGTGAGTCCTCACCACTTATGGCTCCATAGCCTGCTATTTAAGTGGGTTACCGGTCTCCGCCAAGTAGCTGGTGTAAGAACACAGTAACTGAGCCCAGTGTGATCAGCCAAAAAGAGGTACCAATCATCATAAACCATCCCTTGAGTCTCGGTTCTGCTGGTTTCGGACGTCGTGTGGATGCAGCGGATTTGACTACCGTCCTCATAGGAATGCCGGATGTAATGAAACTTCCGCTTCCAAATATACGATATCAAAGGTTGTGGAGACGATGGGCGAACTTGGCAGCGGCCCCCCCCACCGGGGGTCTCCGGCGTAGCGGTACGGTCTATGCTAAGGCGTTGCTAACAATTGCAGGAGCACGGGGCTCGCAAGTAAAAGCACCGCACTGGGCATGATACCGGGGAATACGGAGTCTTCCCTTATGCCGAAAGAAGCAGCTATAACTTCCTCGGGTAAAGGGCAAGAGAAGATCGTAGGCACGTACTCCCGAACTTCAAGAGATCCCGGTTTGCTGCGCCAACCCAATGGTAGCCACATCACGCATATTAGACCGTTGCTGAAATAGTAAAGGCCGCAACCTTCAGATGTCAGCCTTTTCATGCTGTGGATTAACAAGAGTGGGGAAGCAATACGAAGTGAGTTCGTTGGGCATGCGGGAGGGCGGCAGGAAGCAAACGGGTTGCGGCCCGGCGCGGTACGTTGTGAATCGATCTCTGACGCATACCCTCCAGCAATTCCTAAAACCTCCGCATTTTTAGATGTCTGCTGTCGGTCAGGTAGTCAACAGGTTTGTTCACCGAAACGACTGGTCTTCACCCCGTCAAATCATTAAACGCGCCCGCAGTGCTTTCACGGGTCCCCGACGTCAGATCGCCCTAGACCATGATGCCCGGTACCAAAGTCTCACTGCCGTCACGGTAAGTGGTATATGCGGTTGGGCGGCTCTCTACTTCGGTTGATGAATAATGGTGCTGAAGGCGACT";
  haplotypes[4].align_begin_wrt_ref = 0;
  haplotypes[4].score = -1.63548;
  haplotypes[4].cigar = "299M1I725M";

  haplotypes[5].seq = "AATCGAAGGTCGTAAGGACACGGTTGAGCGTTCAGCGTTCATGTGAGTCCTCACCACTTATGGCTCCATAGCCTGCTATTTAAGTGGGTTACCGGTCTCCGCCAAGTAGCTGGTGTAAGAACACAGTAACTGAGCCCAGTGTGATCAGCCTTTTTGAGGTACCAATCATCATAAACCATCCCTTGAGTCTCGGTTCTGCTGGTTTCGGACGTCGTGTGGATGCAGCGGATTTGACTACCGTCCTCATAGGAATGCCGGATGTAATGAAACTTCCGCTTCCAAATATACGATATCAAAGGTGTGGAGACGATGGGCGAACTTGGCAGCGGCCCCCCCCACCGGGGGTCTCCGGCGTAGCGGTACGGTCTATGCTAAGGCGTTGCTAACAATTGCAGGAGCACGGGGCTCGCAAGTAAAAGCACCGCACTGGGCATGATACCGGGGAATACGGAGTCTTCCCTTATGCCGAAAGAAGCAGCTATAACTTCCTCGGGTAAAGGGCAAGAGAAGATCGTAGGCACGTACTCCCGAACTTCAAGAGATCCCGGTTTGCTGCGCCAACCCAATGGTAGCCACATCACGCATATTAGACCGTTGCTGAAATAGTAAAGGCCGCAACCTTCAGATGTCAGCCTTTTCATGCTGTGGATTAACAAGAGTGGGGAAGCAATACGAAGTGAGTTCGTTGGGCATGCGGGAGGGCGGCAGGAAGCAAACGGGTTGCGGCCCGGCGCGGTACGTTGTGAATCGATCTCTGACGCATACCCTCCAGCAATTCCTAAAACCTCCGCATTTTTAGATGTCTGCTGTCGGTCAGGTAGTCAACAGGTTTGTTCACCGAAACGACTGGTCTTCACCCCGTCAAATCATTAAACGCGCCCGCAGTGCTTTCACGGGTCCCCGACGTCAGATCGCCCTAGACCATGATGCCCGGTACCAAAGTCTCACTGCCGTCACGGTAAGTGGTATATGCGGTTGGGCGGCTCTCTACTTCGGTTGATGAATAATGGTGCTGAAGGCGACT";
  haplotypes[5].align_begin_wrt_ref = 0;
  haplotypes[5].score = -1.63548;
  haplotypes[5].cigar = "150M5D5I869M";

  haplotypes[6].seq = "AATCGAAGGTCGTAAGGACACGGTTGAGCGTTCAGCGTTCATGTGAGTCCTCACCACTTATGGCTCCATAGCCTGCTATTTAAGTGGGTTACCGGTCTCCGCCAAGTAGCTGGTGTAAGAACACAGTAACTGAGCCCAGTGTGATCAGCCAAAAAGAGGTACCAATCATCATAAACCATCCCTTGAGTCTCGGTTCTGCTGGTTTCGGACGTCGTGTGGATGCAGCGGATTTGACTACCGTCCTCATAGGAATGCCGGATGTAATGAAACTTCCGCTTCCAAATATACGATATCAAAGGTGTGGAGACGATGGGCGAACTTGGCAGCGGCCCCCCCCACCGGGGGTCTCCGGCGTAGCGGTACGGTCTATGCTAAGGCGTTGCTAACAATTGCAGGAGCACGGGGCTCGCAAGTAAAAGCACCGCACTGGGCATGATACCGGGGAATACGGAGTCTTCCCTTATGCCGAAAGAAGCAGCTATAACTTCCTCGGGTAAAGGGCAAGAGAAGATCGTAGGCACGTACTCCCGAACTTCAAGAGATCCCGGTTTGCTGCGCCAACCCAATGGTAGCCACATCACGCATATTAGACCGTTGCTGAAATAGTAAAGGCCGCAACCTTCAGATGTCAGCCTTTTCATGCTGTGGATTAACAAGAGTGGGGAAGCAATACGAAGTGAGTTCGTTGGGCATGCGGGAGGGCGGCAGGAAGCAAACGGGTTGCGGCCCGGCGCGGTACGTTGTGAATCGATCTCTGACGCATACCCTCCAGCAATTCCTAAAACCTCCGCATTTTTAGTGTCTGCTGTCGGTCAGGTAGTCAACAGGTTTGTTCACCGAAACGACTGGTCTTCACCCCGTCAAATCATTAAACGCGCCCGCAGTGCTTTCACGGGTCCCCGACGTCAGATCGCCCTAGACCATGATGCCCGGTACCAAAGTCTCACTGCCGTCACGGTAAGTGGTATATGCGGTTGGGCGGCTCTCTACTTCGGTTGATGAATAATGGTGCTGAAGGCGACT";
  haplotypes[6].align_begin_wrt_ref = 0;
  haplotypes[6].score = -1.63548;
  haplotypes[6].cigar = "799M1D224M";

  haplotypes[7].seq = "AATCGAAGGTCGTAAGGACACGGTTGAGCGTTCAGCGTTCATGTGAGTCCTCACCACTTATGGCTCCATAGCCTGCTATTTAAGTGGGTTACCGGTCTCCGCCAAGTAGCTGGTGTAAGAACACAGTAACTGAGCCCAGTGTGATCAGCCAAAAAGAGGTACCAATCATCATAAACCATCCCTTGAGTCTCGGTTCTGCTGGTTTCGGACGTCGTGTGGATGCAGCGGATTTGACTACCGTCCTCATAGGAATGCCGGATGTAATGAAACTTCCGCTTCCAAATATACGATATCAAAGGTGTGGAGACGATGGGCGAACTTGGCAGCGGCCCCCCCCACCGGGGGTCTCCGGCGTAGCGGTACGGTCTATGCTAAGGCGTTGCTAACAATTGCAGGAGCACGGGGCTCGCAAGTAAAAGCACCGCACTGGGCATGATACCGGGGAATACGGAGTCTTCCCTTATGCCGAAAGAAGCAGCTATAACTTCCTCGGGTAAAGGGCAAGAGAAGATCGTAGGCACGTACTCCCGAACTTCAAGAGATCCCGGTTTGCTGCGCCAACCCAATGGTAGCCACATCACGCATATTAGACCGTTGCTGAAATAGTAAAGGCCGCAACCTTCAGATGTCAGCCTTTTCATGCTGTGGATTAACAAGAGTGGGGAAGCAATACGAAGTGAGTTCGTTGGGCATGCGGGAGGGCGGCAGGAAGCAAACGGGTTGCGGCCCGGCGCGGTACGTTGTGAATCGATCTCTGACGCATACCCTCCAGCAATTCCTAAAACCTCCGCATTTTTAGATGTCTGCTGTCGGTCAGGTAGTCAACAGGTTTGTTCACCGAAACGACTGGTCTTCACCCCGTCAAATCATTAAACGCGCCCGCAGTGCTTTCACGGGTCCCCGACGTCAGATCGCCCTAGACCATGATGCCCGGTACCAAAGTCTCACTGCCGTCACGGTAAGTGGTATATGCGGTTGGGCGGCTCTCTACTTCGGTTGATGAATAATGGTGCTGAAGGCGACT";
  haplotypes[7].align_begin_wrt_ref = 0;
  haplotypes[7].score = -2.33445;
  haplotypes[7].cigar = "1024M";

  auto pairhmm = PairHMM{};
  auto likelihoods = pairhmm.compute_likelihoods(haplotypes, reads);

  std::vector<std::vector<double>> ans = {
    {-3.010440, -7.510021, -3.010021, -3.010864, -7.510021, -3.010440, -7.510021, -7.510021},
    {-3.010440, -7.510021, -3.010021, -3.010864, -7.510021, -3.010440, -7.510021, -7.510021},
    {-3.010440, -7.510021, -3.010021, -3.010864, -7.510021, -3.010440, -7.510021, -7.510021},
    {-3.010440, -7.510021, -3.010021, -3.010864, -7.510021, -3.010440, -7.510021, -7.510021},
    {-3.010440, -7.510021, -3.010021, -3.010864, -7.510021, -3.010440, -7.510021, -7.5100212}};

  for (auto i = 0; i < likelihoods.size(); i++)
    for (auto j = 0; j < likelihoods[i].size(); j++)
      REQUIRE(likelihoods[i][j] == Approx(ans[i][j]));
}
