#include <biovoltron/file_io/reference.hpp>
#include <iostream> //debug
#include <catch.hpp>
#include <filesystem>
#include <sstream>


TEST_CASE("Reference") {
  SECTION("origin_seq") {
    biovoltron::ReferenceRecord<false> ref_record;
    biovoltron::ReferenceRecord<true> ref_record_encoded;
    SECTION("case: with unknown intervals") {
      ref_record.seq = "ACTGAACAGTCCATCGTGACTGGACTGACTGA";
      ref_record.unknown_intervals = {{5, 7}, {12, 15}};
      auto o_seq = ref_record.origin_seq();
      REQUIRE(o_seq == "ACTGANNAGTCCNNNGTGACTGGACTGACTGA");
    }
    SECTION("case: without unknown intervals") {
      ref_record.seq = "ACTGAACAGTCCATCGTGACTGGACTGACTGA";
      ref_record.unknown_intervals = {};
      auto o_seq = ref_record.origin_seq();
      REQUIRE(o_seq == "ACTGAACAGTCCATCGTGACTGGACTGACTGA");
    }
    SECTION("case: with unknown intervals at begin and end") {
      ref_record.seq = "ACTGAACAGTCCATCGTGACTGGACTGACTGA";
      ref_record.unknown_intervals = {{0, 2}, {30, 32}};
      auto o_seq = ref_record.origin_seq();
      REQUIRE(o_seq == "NNTGAACAGTCCATCGTGACTGGACTGACTNN");
    }
    SECTION("case: with unknown intervals, Encoded = true") {
      ref_record_encoded.seq = biovoltron::istring{0, 1, 2, 3, 0, 0, 1, 2, 3, 3, 1, 2};
      ref_record_encoded.unknown_intervals = {{2,5}, {8, 10}};
      auto o_seq = ref_record_encoded.origin_seq();
      auto expected = biovoltron::istring{0, 1, 4, 4, 4, 0, 1, 2, 4, 4, 1, 2};
      REQUIRE(o_seq == expected);
    }
    SECTION("case: without unknown intervals, Encoded = true") {
      ref_record_encoded.seq = biovoltron::istring{0, 1, 2, 3, 0, 0, 1, 2, 3, 3, 1, 2};
      ref_record_encoded.unknown_intervals = {};
      auto o_seq = ref_record_encoded.origin_seq();
      auto expected = biovoltron::istring{0, 1, 2, 3, 0, 0, 1, 2, 3, 3, 1, 2};
      REQUIRE(o_seq == expected);
    }
    SECTION("case: with unknown intervals at begin and end, Encoded = true") {
      ref_record_encoded.seq = biovoltron::istring{0, 1, 2, 3, 0, 0, 1, 2, 3, 3, 1, 2};
      ref_record_encoded.unknown_intervals = {{0,2}, {10,12}};
      auto o_seq = ref_record_encoded.origin_seq();
      auto expected = biovoltron::istring{4, 4, 2, 3, 0, 0, 1, 2, 3, 3, 4, 4};
      REQUIRE(o_seq == expected);
    }
    SECTION("case: with invalid unknown_intervals") {
      ref_record.seq = "ACGT";
      ref_record.unknown_intervals = {{2, 10}};
      REQUIRE_THROWS_AS(ref_record.origin_seq(), std::out_of_range);
    }

  }
  SECTION("save and load") {
    biovoltron::ReferenceRecord<false> ref_record;
    biovoltron::ReferenceRecord<true> ref_record_encoded;
    SECTION("case: Encoded = false -> Encoded = false") {
      const std::string filename = "test_ref_false.bfa";
      auto ref = biovoltron::ReferenceRecord<false>{};
      ref.species = "TestSpecies";
      ref.chr_num = 2;
      ref.chr_names = {"chr1", "chr2"};
      ref.seq = "ACGTNNACGT";
      ref.base_cnt = {2, 2, 2, 2, 2}; // A:2 C:2 G:2 T:2 N:2
      ref.chr_end_pos = {5, 10};
      ref.unknown_intervals = {{4, 6}}; // N at 4,5

      // Save to file
      {
        auto fout = std::ofstream{filename, std::ios::binary};
        REQUIRE(fout.is_open());
        ref.save(fout);
        fout.close();
      }
      // Load from file
      auto loaded = biovoltron::ReferenceRecord<false>{};
      {
        std::ifstream fin(filename, std::ios::binary);
        REQUIRE(fin.is_open());
        loaded.load(fin);
        fin.close();
      }

      REQUIRE(ref.species == loaded.species);
      REQUIRE(ref.chr_num == loaded.chr_num);
      REQUIRE(ref.chr_names == loaded.chr_names);
      REQUIRE(ref.origin_seq() == loaded.origin_seq());
      REQUIRE(ref.base_cnt == loaded.base_cnt);
      REQUIRE(ref.chr_end_pos == loaded.chr_end_pos);
      REQUIRE(ref.unknown_intervals == loaded.unknown_intervals);
      
      REQUIRE(std::filesystem::remove(filename));

    }
    SECTION("case: Encoded = true -> Encoded = true") {
      const std::string filename_encoded = "test_ref_true.bfa";

      ref_record_encoded.species = "TestSpecies";
      ref_record_encoded.chr_num = 2;
      ref_record_encoded.chr_names = {"chr1", "chr2"};
      ref_record_encoded.seq = biovoltron::istring{0, 1, 2, 3, 4, 4, 0, 1, 2, 3};
      ref_record_encoded.base_cnt = {2, 2, 2, 2, 2}; // A:2 C:2 G:2 T:2 N:2
      ref_record_encoded.chr_end_pos = {5, 10};
      ref_record_encoded.unknown_intervals = {{4, 6}}; // N at 4,5

      // Save to file
      {
        std::ofstream fout(filename_encoded, std::ios::binary);
        REQUIRE(fout.is_open());
        ref_record_encoded.save(fout);
        fout.close();
      }

      // Load from file
      biovoltron::ReferenceRecord<true> loaded;
      {
        std::ifstream fin(filename_encoded, std::ios::binary);
        REQUIRE(fin.is_open());
        loaded.load(fin);
        fin.close();
      }
      
      REQUIRE(ref_record_encoded.species == loaded.species);
      REQUIRE(ref_record_encoded.chr_num == loaded.chr_num);
      REQUIRE(ref_record_encoded.chr_names == loaded.chr_names);
      REQUIRE(ref_record_encoded.origin_seq() == loaded.origin_seq());
      REQUIRE(ref_record_encoded.base_cnt == loaded.base_cnt);
      REQUIRE(ref_record_encoded.chr_end_pos == loaded.chr_end_pos);
      REQUIRE(ref_record_encoded.unknown_intervals == loaded.unknown_intervals);

      REQUIRE(std::filesystem::remove(filename_encoded));
    }
    SECTION("case: Encoded = false -> Encoded = true") {
      const std::string filename_mixed = "test_ref_mixed.bfa";
      ref_record.species = "TestSpecies";
      ref_record.chr_num = 2;
      ref_record.chr_names = {"chr1", "chr2"};
      ref_record.seq = "ACGTNNACGT";
      ref_record.base_cnt = {2, 2, 2, 2, 2}; // A:2 C:2 G:2 T:2 N:2
      ref_record.chr_end_pos = {5, 10};
      ref_record.unknown_intervals = {{4, 6}}; // N at 4,5

      // Save to file
      {
        std::ofstream fout(filename_mixed, std::ios::binary);
        REQUIRE(fout.is_open());
        ref_record.save(fout);
        fout.close();
      }

      // Load into encoded record
      biovoltron::ReferenceRecord<true> loaded_encoded;
      {
        std::ifstream fin(filename_mixed, std::ios::binary);
        REQUIRE(fin.is_open());
        loaded_encoded.load(fin);
        fin.close();
      }
      
      REQUIRE(ref_record.species == loaded_encoded.species);
      REQUIRE(ref_record.chr_num == loaded_encoded.chr_num);
      REQUIRE(ref_record.chr_names == loaded_encoded.chr_names);
      REQUIRE(loaded_encoded.origin_seq() == biovoltron::istring{0, 1, 2, 3, 4, 4, 0, 1, 2, 3});
      REQUIRE(ref_record.base_cnt == loaded_encoded.base_cnt);
      REQUIRE(ref_record.chr_end_pos == loaded_encoded.chr_end_pos);
      REQUIRE(ref_record.unknown_intervals == loaded_encoded.unknown_intervals);

      REQUIRE(std::filesystem::remove(filename_mixed));
    }
    SECTION("case: Encoded = true -> Encoded = false") {
      const std::string filename_mixed_encoded = "test_ref_mixed_encoded.bfa";
      ref_record_encoded.species = "TestSpecies";
      ref_record_encoded.chr_num = 2;
      ref_record_encoded.chr_names = {"chr1", "chr2"};
      ref_record_encoded.seq = biovoltron::istring{0, 1, 2, 3, 4, 4, 0, 1, 2, 3};
      ref_record_encoded.base_cnt = {2, 2, 2, 2, 2}; // A:2 C:2 G:2 T:2 N:2
      ref_record_encoded.chr_end_pos = {5, 10};
      ref_record_encoded.unknown_intervals = {{4, 6}}; // N at 4,5  
      // Save to file
      {
        std::ofstream fout(filename_mixed_encoded, std::ios::binary);
        REQUIRE(fout.is_open());
        ref_record_encoded.save(fout);
        fout.close();
      }
      // Load into non-encoded record
      biovoltron::ReferenceRecord<false> loaded_non_encoded;
      {
        std::ifstream fin(filename_mixed_encoded, std::ios::binary);
        REQUIRE(fin.is_open());
        loaded_non_encoded.load(fin);
        fin.close();  
      }
      REQUIRE(ref_record_encoded.species == loaded_non_encoded.species);
      REQUIRE(ref_record_encoded.chr_num == loaded_non_encoded.chr_num);
      REQUIRE(ref_record_encoded.chr_names == loaded_non_encoded.chr_names);
      REQUIRE(loaded_non_encoded.origin_seq() == "ACGTNNACGT");
      REQUIRE(ref_record_encoded.base_cnt == loaded_non_encoded.base_cnt);
      REQUIRE(ref_record_encoded.chr_end_pos == loaded_non_encoded.chr_end_pos);
      REQUIRE(ref_record_encoded.unknown_intervals == loaded_non_encoded.unknown_intervals);
      REQUIRE(std::filesystem::remove(filename_mixed_encoded));
  
    }
  }
  SECTION("operator>>") {
    SECTION("case: Encoded = false") {
      const std::string input_data = ">chr1\n"
                                    "ACGTNNAC\n"
                                    ">chr2\n"
                                    "GANN\n";
      std::istringstream iss(input_data);
      biovoltron::ReferenceRecord<false> ref_record;
      iss >> ref_record;
      REQUIRE(ref_record.chr_num == 2);
      REQUIRE(ref_record.chr_names == std::vector<std::string>{"chr1", "chr2"});
      REQUIRE(ref_record.origin_seq() == "ACGTNNACGANN");
      REQUIRE(ref_record.chr_end_pos == std::vector<std::uint32_t>{8, 12});
      REQUIRE(ref_record.base_cnt == std::vector<std::uint32_t>{3, 2, 2, 1, 4}); // A,C,G,T,N counts
      REQUIRE(ref_record.unknown_intervals == std::vector<std::array<std::uint32_t, 2>>{
        {4, 6}, // NN in chr1
        {10, 12} // NN in chr2
      });
    }
    SECTION("case: Encoded = true") {
      const std::string fasta = 
        ">X\n"
        "ACNNGT";

      std::istringstream iss(fasta);
      auto record = biovoltron::ReferenceRecord<true>{};
      iss >> record;

      REQUIRE(record.chr_num == 1);
      REQUIRE(record.chr_names == std::vector<std::string>{"X"});
      REQUIRE(record.origin_seq() == biovoltron::istring{0, 1, 4, 4, 2, 3});
      REQUIRE(record.chr_end_pos == std::vector<std::uint32_t>{6});
      REQUIRE(record.base_cnt == std::vector<std::uint32_t>{1, 1, 1, 1, 2}); // A,C,G,T,N
      REQUIRE(record.unknown_intervals == std::vector<std::array<std::uint32_t, 2>>{
        {2, 4}
      });
    }
    SECTION("case: with long sequence") {
      std::string input = ">chr1\n" + std::string(1000, 'A') +
                          "\n>chr2\n" + std::string(1000, 'C');
      std::istringstream iss(input);
      biovoltron::ReferenceRecord<false> ref;
      iss >> ref;
      
      REQUIRE(ref.chr_num == 2);
      REQUIRE(ref.chr_names == std::vector<std::string>{"chr1", "chr2"});
      REQUIRE(ref.origin_seq() == std::string(1000, 'A') + std::string(1000, 'C'));
      REQUIRE(ref.chr_end_pos == std::vector<uint32_t>{1000, 2000});
      REQUIRE(ref.base_cnt == std::vector<std::uint32_t>{1000, 1000, 0, 0, 0});
      REQUIRE(ref.unknown_intervals.empty());
    }
  }
}


