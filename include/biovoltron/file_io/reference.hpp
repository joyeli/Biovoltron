#pragma once

#include <biovoltron/utility/istring.hpp>
#include <biovoltron/file_io/fasta.hpp>
#include <biovoltron/utility/archive/serializer.hpp>
#include <biovoltron/utility/istring.hpp>
#include <experimental/random>

namespace biovoltron {


/*
 * A struct that store reference of a species
 * which implemented using C++20
 *
 * sequence could be store as two types, DibitVector<std::uint8_t>, or std::string
 * depend on the variable Encoded
 *
 * Example
 * ```cpp
 * #include <iostream>
 * #include <experimental/random>
 * 
 * using namespace biovoltron;
 *
 * int main(void) {
 *   auto ref = ReferenceRecord<false> { .species = "Human" };
 *   std::ifstream {"GRCh38"} >> ref;
 *
 *   //chromosome name
 *   for(const auto &name : ref.chr_names) {
 *     std::cout << name << "\n";
 *   }
 *
 *   //base count
 *   for(const auto &cnt : ref.base_cnt) {
 *     std::cout << cnt << " ";
 *   }
 *
 *   //end position in each chromosome
 *   for(const auto &pos : ref.chr_end_pos) {
 *     std::cout << pos << " ";
 *   }
 *
 *   //unknown intervals
 *   for(const auto &row : ref.unknown_intervals) [
 *     std::cout << row[0] << " - " << row[1] << "\n";
 *   }
 *
 *   //save
 *   {
 *     auto fout = std::ostream{"GRCh38.bfa", std::ios::binary};
 *     ref.save(fout);
 *     fout.close();
 *   }
 *
 *   //load
 *   auto load_ref = ReferenceRecord<false>{};
 *   {
 *     auto fin = std::ifstream{"GRCh38.bfa", std::ios::binary};
 *     load.ref.load(fin);
 *     fin.close();
 *   }
 *
 *   assert(ref == load_ref);
 *
 * }
 * ```
 *
 */


template<bool Encoded = false>
struct ReferenceRecord {
  constexpr static auto encoded = Encoded;
  std::string species {};
  std::int32_t chr_num {};

  std::vector<std::string> chr_names {};
  std::conditional_t<Encoded, istring, std::string> seq {};
  
  char (*substitute)(char&) = [](char &ch) { return (ch == 'N' or ch == 'n') ? 
      Codec::to_char(std::experimental::randint(0, 3)) : static_cast<char>(std::toupper(ch)); };

  std::int8_t (*substitute_encoded)(std::int8_t&) = [](std::int8_t &ch) { return (ch > 3) ? 
      static_cast<std::int8_t>(std::experimental::randint(0, 3)) : ch; };
  
  std::vector<std::uint32_t> base_cnt{}; //0:A, 1:C, 2:G, 3:T, 4:N
  std::vector<std::uint32_t> chr_end_pos {};
  std::vector<std::array<std::uint32_t, 2>> unknown_intervals {}; //[begin, end)

  bool
  operator==(const ReferenceRecord &other) const = default;

  auto
  origin_seq() {
    std::conditional_t<Encoded, istring, std::string> o_seq {};
    auto start = std::uint32_t {}, len = std::uint32_t{} ;
    for( auto n_it = unknown_intervals.begin(); 
         n_it != unknown_intervals.end(); ++n_it ) {
      if( start <= (*n_it)[0] ) {
        o_seq += seq.substr(start, (*n_it)[0] - start );
      }
      {
        if constexpr ( encoded ) {
          o_seq += Codec::to_istring(std::string( (*n_it)[1] - (*n_it)[0], 'N' ));
        }
        else {
          o_seq += std::string( (*n_it)[1] - (*n_it)[0], 'N' );
        }
      }
      start = (*n_it)[1];
    }
    o_seq += seq.substr(start);
    return o_seq;
  }

  auto
  save(std::ofstream &fout) const {
    auto Dibit = DibitVector<std::uint8_t> {};
    Dibit.reserve( seq.size() );
    for( const auto &ch : seq ) {
      if constexpr (encoded){
        Dibit.push_back( ch );
      }
      else {
        Dibit.push_back( Codec::to_int(ch) );
      }
    }
    fout.write(reinterpret_cast<const char*>(&chr_num), sizeof(chr_num));
    Serializer::save(fout, species);
    Serializer::save(fout, Dibit);
    Serializer::save(fout, base_cnt);
    Serializer::save(fout, chr_end_pos);
    Serializer::save(fout, unknown_intervals);

    for(auto i {0}; i < chr_num; ++i) {
      Serializer::save(fout, chr_names[i]);
    }
  }

  auto
  load(std::ifstream &fin) {
    auto Dibit = DibitVector<std::uint8_t> {};
    fin.read(reinterpret_cast<char*>(&chr_num), sizeof(chr_num));
    Serializer::load(fin, species);
    Serializer::load(fin, Dibit);
    Serializer::load(fin, base_cnt);
    Serializer::load(fin, chr_end_pos);
    Serializer::load(fin, unknown_intervals);
    chr_names = std::vector<std::string> (chr_num, "");
    for(auto i {0}; i < chr_num; ++i) {
      auto name = std::string {};
      Serializer::load(fin, name);
      chr_names[i] = name;
    }
    seq.reserve( Dibit.size() );
    auto base_view = std::views::transform([](auto c){ if constexpr (encoded) return c; else  return "ACGT"[c]; });
    for( const auto &ch : Dibit | base_view ){
      seq.push_back( ch );
    }
  }

};

template<bool Encoded>
inline auto&
operator>>(std::istream &is, ReferenceRecord<Encoded> &record) {
  record.base_cnt = std::vector<std::uint32_t> (5, 0);
  auto ref = FastaRecord<Encoded>{};
  auto last_chr_pos = std::uint32_t {0};
  std::conditional_t<Encoded, std::int8_t, char> ch {};
  while (is >> ref) {
    ++record.chr_num;
    record.chr_end_pos.push_back( last_chr_pos + ref.seq.size() );
    record.chr_names.push_back( ref.name );
    for( auto i = std::uint32_t {0}; i < ref.seq.size(); ++i ) {
      if( ref.seq[i] == 'N' or ref.seq[i] == 'n' or ref.seq[i] == 4) {
        if( i == 0 or ref.seq[i] != ch ) {
          record.unknown_intervals.push_back( {last_chr_pos + i, last_chr_pos + i + 1} );
        }
        else {
          ++record.unknown_intervals.back()[1];
        }
        ++record.base_cnt[4];
      }
      else { // not n
        if constexpr (Encoded) {
          ++record.base_cnt[ref.seq[i]];
        }
        else {
          ++record.base_cnt[Codec::to_int(ref.seq[i])];
        }
      }
      ch = ref.seq[i];
    }
    last_chr_pos = record.chr_end_pos.back();
    if constexpr (Encoded) {
      std::transform(ref.seq.begin(), ref.seq.end(), ref.seq.begin(), record.substitute_encoded );
    }
    else {
      std::transform(ref.seq.begin(), ref.seq.end(), ref.seq.begin(), record.substitute );
    }
    //std::transform(ref.seq.begin(), ref.seq.end(), ref.seq.begin(),
    //  [](unsigned char c) { return std::toupper(c); } );

    for(const auto &ch : ref.seq) {
      record.seq.push_back(ch);
    }
      
  }
  return is;
}
}  // namespace biovoltron
