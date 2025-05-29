#pragma once

#include <htslib/sam.h>
#include <spdlog/spdlog.h>
#include <biovoltron/file_io/sam.hpp>
#include <filesystem>
#include <vector>
#include <cstring>
#include <iostream>
#include <limits>
#include <map>
#include <fstream>

namespace biovoltron {

/**
 * @ingroup file_io
 * @brief A class about the config of SAM file and some core methods
 */

struct IBamStream : protected std::ifstream {

  IBamStream(const std::filesystem::path&);
  ~IBamStream();
  
  friend auto& operator>>(IBamStream& is, SamHeader& h);
  template<bool Encoded>
  friend auto& operator>>(IBamStream& is, SamRecord<Encoded>& r);

  explicit operator bool() const { return !is_eof; };

  auto open(const std::filesystem::path&);
  auto is_open();
  auto eof();
  auto close();
  auto load_index(const std::filesystem::path&);
  auto set_region(std::string_view rname, std::size_t beg, std::size_t end);
  auto is_indexed() const noexcept;
  auto on_sequential() const noexcept;
  auto on_ranged() const noexcept;
  auto on_unmapped() const noexcept;
  auto set_unmapped();
  auto set_sequential();
  auto to_begin();

 private:
  enum class IterateType : std::uint8_t { SEQUENTIAL, RANGED, UNMAPPED };

  auto traverse_on(std::int32_t tid, std::size_t beg, std::size_t end, IterateType type);
  auto clear();

  samFile* bam = nullptr;
  bam_hdr_t* bam_header = nullptr;
  hts_idx_t* idx = nullptr;
  hts_itr_t* itr = nullptr;
  bool is_eof = false;
  std::filesystem::path path;
  IterateType itr_type = IterateType::SEQUENTIAL;
};

struct OBamStream : protected std::ofstream {
  explicit OBamStream(const std::filesystem::path& path, bool gen_idx = false);
  ~OBamStream();

  friend auto& operator<<(OBamStream& os, SamHeader& h);
  template<bool Encoded>
  friend auto& operator<<(OBamStream& os, SamRecord<Encoded>& r);

  auto open(const std::filesystem::path& bam_path, bool gen_idx);
  auto is_open();
  auto close();

 private:
  auto clear();
  template <typename T, typename U>
  static auto convert_type_and_copy(U source, std::uint8_t* dest);
  template<bool Encoded>
  static auto aux_to_bin(bam1_t* aln, SamRecord<Encoded>& r);
  template<bool Encoded>
  static auto bam_set1(
    bam1_t* aln, SamRecord<Encoded>& r, int tid, std::vector<std::uint32_t>& cigars,
    std::vector<char>& quals, int mtid);

  samFile* bam = nullptr;
  bam_hdr_t* bam_header = nullptr;
  std::map<std::string, std::int32_t> ref_table;
  bool write_idx = false;
  std::filesystem::path path;
  std::filesystem::path idx_path;
  static const inline auto bam_op_table =
    std::map<char, std::uint8_t>{{'M', 0}, {'I', 1}, {'D', 2}, {'N', 3}, {'S', 4},
                                 {'H', 5}, {'P', 6}, {'=', 7}, {'X', 8}, {'B', 9}};
};



/* ----------- IBamStream ----------- */

auto IBamStream::clear() {
  sam_itr_destroy(itr);
  hts_idx_destroy(idx);
  bam_hdr_destroy(bam_header);
  if (bam != nullptr) {
    sam_close(bam);
  }
  bam = nullptr;
  bam_header = nullptr;
  idx = nullptr;
  itr = nullptr;
  is_eof = false;
  itr_type = IterateType::SEQUENTIAL;
}

/**
 * @brief opens a BAM file and associate it with this object. If index file is
 * alongside with the BAM file, it will be open automatically. You can also open
 * the specific index file with `load_index` function
 * 
 * @param bam_path the path to the file to be opened
 */
auto IBamStream::open(const std::filesystem::path& bam_path) {
  if (bam) clear();
  path = bam_path;
  bam = sam_open(bam_path.c_str(), "rb");
  bam_header = sam_hdr_read(bam);
  if (bam_header == nullptr) {
    SPDLOG_WARN("Fail to open header at {}", bam_path.string());
  }
  // load the index
  idx = sam_index_load(bam, bam_path.c_str());
  itr = sam_itr_queryi(idx, HTS_IDX_REST, 0, 0);
}


/**
 * @brief gets the index file handler
 *
 * @return the index file handler, nullptr if there's no index file
 */
auto IBamStream::is_indexed() const noexcept { return idx; };

/**
 * @brief Load or stream a BAM index file
 *
 * @param index_path the path of index file
 * @return whether the index file is been opened correctly
 */
auto IBamStream::load_index(const std::filesystem::path& index_path) {
  idx = sam_index_load2(bam, path.c_str(), index_path.c_str());
  return !!(is_indexed());
}

/**
 * @brief check whether the file is opened
 * 
 * @return true if the file is opened, false otherwise
 */
auto IBamStream::is_open() {
  return this->bam != nullptr;
}

/**
 * @brief check whether the eof(end of file) flag is set
 * 
 * @return true if the eof flag is set, false otherwise
 */
auto IBamStream::eof() {
  return this->is_eof;
}

/**
 * @brief Close the file stream by clear the handler of htslib.
 * 
 * @return `void`
 */
auto IBamStream::close() {
  clear();
}


/**
 * @brief constructs the IBamStream object
 * 
 * @param path the path of the file to be opened
 */
IBamStream::IBamStream(const std::filesystem::path& bam_path) { 
  open(bam_path); 
};

/**
 * @brief destruct the IBamStream object
 */
IBamStream::~IBamStream() { clear(); }


/**
 * @brief checks whether the iterate type is sequential
 * 
 * @return true if the iterator type is SEQUENTIAL, false otherwise
 */
auto IBamStream::on_sequential() const noexcept { return itr_type == IterateType::SEQUENTIAL; };

/**
 * @brief checks whether the iterate type is ranged
 * 
 * @return true if the iterator type is RANGED, false otherwise
 */
auto IBamStream::on_ranged() const noexcept { return itr_type == IterateType::RANGED; };

/**
 * @brief checks whether the iterate type is unmapped
 * 
 * @return true if the iterator type is UNMAPPED, false otherwise
 */
auto IBamStream::on_unmapped() const noexcept { return itr_type == IterateType::UNMAPPED; };

/**
 * @brief setup 
 * 
 * @param tid the tid value for iterate over bam file
 * @param beg the start position
 * @param end the end position
 * @param type the type for iterating file
 * @return auto 
 */
auto IBamStream::traverse_on(std::int32_t tid, std::size_t beg, std::size_t end, IterateType type) {
  if (!is_indexed()) return false;
  auto new_itr = sam_itr_queryi(idx, tid, beg, end);
  if (!new_itr) return false;
  sam_itr_destroy(itr);
  itr = new_itr;
  itr_type = type;
  is_eof = false;
  return true;
}

/**
 * @brief sets the iterator to type RANGED from beg to end
 * 
 * @param rname reference name
 * @param beg start position in target
 * @param end end position in target
 * @return true on success, false otherwise
 */
auto IBamStream::set_region(std::string_view rname, std::size_t beg, std::size_t end) {
  if (int tid = bam_name2id(bam_header, rname.data()); tid >= 0)
    return traverse_on(tid, beg, end, IterateType::RANGED);
  return false;
}

/**
 * @brief set the iterator type to UNMAPPED
 * 
 * @return true on success, false otherwise
 */
auto IBamStream::set_unmapped() {
  return traverse_on(HTS_IDX_NOCOOR, 0, 0, IterateType::UNMAPPED);
};

/**
 * @brief set the iterator type to SEQUENTIAL
 * 
 * @return true on success, false otherwise
 */
auto IBamStream::set_sequential() {
  return traverse_on(HTS_IDX_REST, 0, 0, IterateType::SEQUENTIAL);
};

/**
 * @brief set the iterator back to begin and type to SEQUENTIAL
 * 
 * @return true on success, false otherwise
 */
auto IBamStream::to_begin() { return traverse_on(HTS_IDX_START, 0, 0, IterateType::SEQUENTIAL); };

/**
 * @brief extracts value of the header 
 * 
 * @param is an IBamStream object to be extracted
 * @param h an output SamHeader object
 */
auto& operator>>(IBamStream& is, SamHeader& h) {
  h.lines.clear();
  auto header_str = static_cast<std::string_view>(sam_hdr_str(is.bam_header));
  std::ranges::for_each(header_str | std::ranges::views::split('\n'), [&h](auto const& view) {
    std::string buf;
    for (const auto element : view) buf.push_back(element);
    h.lines.emplace_back(std::move(buf));
  });
  return is;
}

/**
 * @brief extract the value
 * 
 * @param is an IBamStream object to be extracted
 * @param r a SamRecord<false> object
 */
template<bool Encoded>
auto& operator>>(IBamStream& is, SamRecord<Encoded>& r) {
  auto c = 0;
  auto aln = bam_init1();
  if (is.itr) c = sam_itr_next(is.bam, is.itr, aln);
  else
    c = sam_read1(is.bam, is.bam_header, aln);

  if (c >= 0) {
    r.qname = bam_get_qname(aln);
    r.flag = aln->core.flag;
    r.rname = aln->core.tid == -1 ? "*" : is.bam_header->target_name[aln->core.tid];
    r.pos = aln->core.pos + 1;
    r.mapq = aln->core.qual;
    r.cigar.clear();
    auto cigar = bam_get_cigar(aln);
    for (auto i = 0u; i < aln->core.n_cigar; i++) {
      r.cigar.emplace_back(bam_cigar_oplen(cigar[i]), BAM_CIGAR_STR[bam_cigar_op(cigar[i])]);
    }
    r.rnext = aln->core.mtid == -1              ? "*"
              : aln->core.mtid == aln->core.tid ? "="
                                                : is.bam_header->target_name[aln->core.mtid];
    r.pnext = aln->core.mpos + 1;
    r.tlen = aln->core.isize;
    r.seq.clear();
    r.qual.clear();
    if (auto len = aln->core.l_qseq; len != 0) {
      auto seq = bam_get_seq(aln);
      auto qual = bam_get_qual(aln);
      r.seq.reserve(len);
      for (auto i = 0; i < len; i++) {
        auto c = seq_nt16_str[bam_seqi(seq, i)];
        if constexpr (Encoded) {
          r.seq += Codec::to_int(c);
        } else {
          r.seq += c;
        }
      }
      if (qual[0] == 255) r.qual = "*";
      else {
        r.qual.reserve(static_cast<std::size_t>(len));
        for (auto i = 0; i < len; i++) r.qual += static_cast<char>(33 + qual[i]);
      }
    }
    // optional fields
    auto optional_str = kstring_t(KS_INITIALIZE);
    auto start_ptr = static_cast<std::uint8_t*>(bam_get_aux(aln));
    auto end_ptr = static_cast<std::uint8_t*>(aln->data + aln->l_data);
    r.optionals.clear();
    while (end_ptr - start_ptr >= 4) {
      // kstring_t will work unexpectedly, when type is 'Z' or 'H'. We need to work around this.
      if (start_ptr[2] == 'Z' || start_ptr[2] == 'H') {
        auto s = reinterpret_cast<char*>(start_ptr);
        while (*(start_ptr++) != 0)
          ;
        auto& x = r.optionals.emplace_back(s, 2);
        x.push_back(':');
        x.push_back(s[2]);
        x.push_back(':');
        x.append(s + 3);
      } else {
        if ((start_ptr = const_cast<uint8_t*>(
               sam_format_aux1(start_ptr, start_ptr[2], start_ptr + 3, end_ptr, &optional_str))))
          r.optionals.emplace_back(ks_c_str(&optional_str));
      }
      ks_clear(&optional_str);
    }
    ks_free(&optional_str);
  } else if (c == -1)
    is.is_eof = true;
  bam_destroy1(aln);
  return is;
}

/* ----------- OBamStream ----------- */

auto OBamStream::clear() {
  if (write_idx) {
    if (auto c = sam_idx_save(bam); c != 0) {
      throw std::runtime_error("An error occured when saving BAM index.");
    }
  }
  if (bam_header) {
    bam_hdr_destroy(bam_header);
  }
  if (bam) {
    sam_close(bam);
  }
  bam = nullptr;
  bam_header = nullptr;
  ref_table.clear();
  path.clear();
}

auto OBamStream::close() {
  clear();
}

/**
 * @brief opens a BAM file and associate it with this object
 * 
 * @param bam_path the path to the file to be opened
 * @param gen_idx whether to generate a index file
 */
auto OBamStream::open(const std::filesystem::path& bam_path, bool gen_idx) {
  if (bam) clear();
  path = bam_path;
  // open BAM for writing
  bam = sam_open(bam_path.c_str(), "wb");
  bam_header = sam_hdr_init();
  write_idx = gen_idx;
  if (write_idx) {
    idx_path = path.string() + ".bai";
  }
}

/**
 * @brief constructs the OBamStream object
 * 
 * @param bam_path the path to the file to be opened
 * @param gen_idx whether to generate a index file
 */
OBamStream::OBamStream(const std::filesystem::path& path, bool gen_idx) { open(path, gen_idx); };

/**
 * @brief destructs the OBamStream object
 */
OBamStream::~OBamStream() { clear(); }

/**
 * @brief Whether the bam file for output is opening
 * 
 * @return bool
 */
auto OBamStream::is_open() {
  return this->bam != nullptr;
}

template <typename T, typename U>
auto OBamStream::convert_type_and_copy(U source, std::uint8_t* dest) {
  auto d = static_cast<T>(source);
  std::memcpy(dest, &d, sizeof(d));
  return sizeof(d);
}


template<bool Encoded>
auto OBamStream::aux_to_bin(bam1_t* aln, SamRecord<Encoded>& r) {
  auto& len = aln->l_data;
  auto aux_data = aln->data;
  for (auto& opt : r.optionals) {
    aux_data[len++] = opt[0];
    aux_data[len++] = opt[1];
    switch (opt[3]) {
      case 'A': {
        aux_data[len++] = 'A';
        aux_data[len++] = opt[5];
        break;
      }
      case 'i': {
        if (opt[5] == '-') {
          auto num = std::stoi(opt.substr(5));
          if (num >= std::numeric_limits<std::int8_t>::min()) {
            aux_data[len++] = 'c';
            len += OBamStream::convert_type_and_copy<std::int8_t>(num, &aux_data[len]);
          } else if (num >= std::numeric_limits<std::int16_t>::min()) {
            aux_data[len++] = 's';
            len += OBamStream::convert_type_and_copy<std::int16_t>(num, &aux_data[len]);
          } else {
            aux_data[len++] = 'i';
            len += OBamStream::convert_type_and_copy<std::int32_t>(num, &aux_data[len]);
          }
        } else {
          auto num = static_cast<uint32_t>(std::stoul(opt.substr(5)));
          if (num <= std::numeric_limits<std::uint8_t>::max()) {
            aux_data[len++] = 'C';
            len += OBamStream::convert_type_and_copy<std::uint8_t>(num, &aux_data[len]);
          } else if (num <= std::numeric_limits<std::uint16_t>::max()) {
            aux_data[len++] = 'S';
            len += OBamStream::convert_type_and_copy<std::uint16_t>(num, &aux_data[len]);
          } else {
            aux_data[len++] = 'I';
            len += OBamStream::convert_type_and_copy<std::uint32_t>(num, &aux_data[len]);
          }
        }
        break;
      }
      case 'f': {
        aux_data[len++] = 'f';
        auto num = std::stof(opt.substr(5));
        std::memcpy(&aux_data[len], &num, sizeof(num));
        len += sizeof(num);
        break;
      }
      case 'Z':
      case 'H': {
        aux_data[len++] = opt[3];
        std::ranges::copy_n(opt.substr(5).data(), opt.size() - 4, &aux_data[len]);  // -4 = + 1 - 5
        len += opt.size() - 4;
        break;
      }
      case 'B': {
        aux_data[len++] = 'B';
        // subtype
        aux_data[len++] = opt[5];
        // preserve space for length
        auto count = &aux_data[len];
        len += sizeof(std::uint32_t);
        // numbers
        auto pos = 0ul;
        auto arr_length = 0u;
        // auto num_str = ""sv;
        auto remaining = std::string_view(opt.begin() + 5, opt.end());
        // auto remaining = static_cast<std::string_view>(opt.substr(5));
        remaining = remaining.substr(2);
        while (!remaining.empty()) {
          pos = remaining.find(',');
          // parse number
          auto num_str = remaining.substr(0, pos);
          if (opt[5] != 'f') {
            auto num = std::stoll(std::string(num_str));
            switch (opt[5]) {
              case 'c': len += convert_type_and_copy<std::int8_t>(num, &aux_data[len]); break;
              case 'C': len += convert_type_and_copy<std::uint8_t>(num, &aux_data[len]); break;
              case 's': len += convert_type_and_copy<std::int16_t>(num, &aux_data[len]); break;
              case 'S': len += convert_type_and_copy<std::uint16_t>(num, &aux_data[len]); break;
              case 'i': len += convert_type_and_copy<std::int32_t>(num, &aux_data[len]); break;
              case 'I': len += convert_type_and_copy<std::uint32_t>(num, &aux_data[len]);
            }
          } else {
            auto num = std::stof(std::string(num_str));
            std::memcpy(&aux_data[len], &num, sizeof(num));
            len += sizeof(num);
          }
          // move to next number
          if (pos != std::string::npos) remaining = remaining.substr(pos + 1);
          else
            remaining = remaining.substr(remaining.size());
          ++arr_length;
        }
        std::memcpy(count, &arr_length, sizeof(arr_length));
      }
    }
  }
  aln->l_data = len;
}


template<bool Encoded = false>
auto OBamStream::bam_set1(
  bam1_t* aln, SamRecord<Encoded>& r, int tid, std::vector<std::uint32_t>& cigars, std::vector<char>& quals,
  int mtid) {
  // determine data size
  auto rlen = hts_pos_t{0};
  auto qname_nuls = std::size_t{(4 - r.qname.size() % 4) % 4};
  if (!(r.flag & BAM_FUNMAP)) rlen = bam_cigar2rlen(cigars.size(), cigars.data());
  if (rlen == 0) rlen = 1;
  auto data_len = std::size_t{
    r.qname.size() + qname_nuls + cigars.size() * 4 + (r.seq.size() + 1) / 2 +
    r.seq.size()};  // read_name + digar + seq + qual
  // realloc space
  std::uint8_t* new_data;
  auto new_m_data = data_len + 65535;
  if ((bam_get_mempolicy(aln) & BAM_USER_OWNS_DATA) == 0)
    new_data = reinterpret_cast<uint8_t*>(std::realloc(aln->data, new_m_data));
  else {
    if ((new_data = reinterpret_cast<uint8_t*>(std::malloc(new_m_data))) != nullptr) {
      if (aln->l_data > 0)
        std::ranges::copy(
          aln->data, aln->data + (aln->l_data < aln->m_data ? aln->l_data : aln->m_data), new_data);
      bam_set_mempolicy(aln, bam_get_mempolicy(aln) & (~BAM_USER_OWNS_DATA));
    }
  }
  aln->data = new_data;
  aln->m_data = new_m_data;
  
  // set each field in aln structure
  aln->l_data = data_len;
  aln->core.pos = r.pos == 0 ? 0 : r.pos - 1;
  aln->core.tid = tid;
  aln->core.bin = hts_reg2bin(aln->core.pos, aln->core.pos + rlen, 14, 5);
  aln->core.qual = r.mapq;
  aln->core.l_extranul = qname_nuls - 1;
  aln->core.flag = r.flag;
  aln->core.l_qname = r.qname.size() + qname_nuls;
  aln->core.n_cigar = cigars.size();
  aln->core.l_qseq = r.seq.size();
  aln->core.mtid = mtid;
  aln->core.mpos = r.pnext == 0 ? -1ll : r.pnext - 1;
  aln->core.isize = r.tlen;

  // qname
  auto data_ptr = aln->data;
  std::ranges::copy(r.qname.begin(), r.qname.end(), data_ptr);
  for (auto i = 0; i < qname_nuls; i++) data_ptr[r.qname.size() + i] = '\0';
  data_ptr += r.qname.size() + qname_nuls;

  // cigar
  if (cigars.size() > 0) std::memcpy(data_ptr, cigars.data(), cigars.size() * 4);
  data_ptr += cigars.size() * 4;
  
  // seq
  {
    if (r.seq == "*") {
      aln->core.l_qseq = 0;
    } else {
      for (auto i = 0u; i + 1 < r.seq.size(); i += 2) {
        auto a = r.seq[i];
        auto b = r.seq[i + 1];
        if constexpr (Encoded) {
          a = Codec::to_char(a);
          b = Codec::to_char(b);
        }
        *data_ptr++ = (seq_nt16_table[a] << 4) | seq_nt16_table[b];
      }
      // odd length
      if (r.seq.size() & 1) {
        auto a = r.seq[r.seq.size() - 1];
        if constexpr (Encoded) {
          a = Codec::to_char(a);
        }
        *data_ptr++ = seq_nt16_table[a] << 4;
      }
      // for (; i + 1 < r.seq.size(); i += 2)
      //   *data_ptr++ = (seq_nt16_table[(unsigned char)r.seq[i]] << 4) |
      //                 seq_nt16_table[(unsigned char)r.seq[i + 1]];
      // for (; i < r.seq.size(); i++) *data_ptr++ = seq_nt16_table[(unsigned char)r.seq[i]] << 4;
    }
  }
  // qual
  if (r.qual != "*") {
   std::ranges::copy(quals.begin(), quals.end(), data_ptr);
  } else if (r.seq != "*") { // ? will this affected by Encoded
    std::ranges::fill_n(data_ptr, quals.size(), '\xff');
  }
  // optional fields
  OBamStream::aux_to_bin(aln, r);
}

/** 
 * @brief inserts SAM file header to the object
 * 
 * @param os an OBamStream object to be inserted
 * @param h an input SamHeader object
 */
auto& operator<<(OBamStream& os, SamHeader& h) {
  auto header = os.bam_header;
  for (const auto& line : h.lines) {
    int r = sam_hdr_add_lines(header, line.c_str(), line.size());
  }
  if (auto c_hdr_write = sam_hdr_write(os.bam, header); c_hdr_write < 0) { /* TODO: error handler */
    throw std::runtime_error("An error occured when writing to a BAM file.");
  } else {
    // Build reference table
    for (auto i = 0; i < header->n_targets; ++i) os.ref_table[header->target_name[i]] = i;
    // build idx
    if (os.write_idx) {
      if (auto c_idx_init = sam_idx_init(os.bam, header, 0, os.idx_path.c_str());
          c_idx_init < 0) { /* TODO: error handler */
        
      };
    }
  }
  return os;
}

/** 
 * @brief inserts SAM file to the object
 * 
 * @param os an OBamStream object to be inserted
 * @param r an input SamRecord<false> object
 */
template<bool Encoded>
auto& operator<<(OBamStream& os, SamRecord<Encoded>& r) {
  // prepare some data
  auto tid = r.rname == "*" ? -1 : os.ref_table[r.rname];
  auto cigars = std::vector<std::uint32_t>{};
  cigars.reserve(r.cigar.size());
  for (const auto& c : r.cigar)
    cigars.emplace_back(bam_cigar_gen(c.size, os.bam_op_table.at(c.op)));
  auto quals = std::vector<char>(r.seq.size(), 255);
  if (r.qual != "*") {
    std::ranges::transform(r.qual, quals.begin(), [](auto c) { return c - 33; });
  }
  auto mtid = r.rnext == "*" ? -1 : r.rnext == "=" ? tid : os.ref_table[r.rnext];
  /*
  bam_set1(
    os.aln,
    r.qname.size(),
    r.qname.data(),
    r.flag,
    tid,
    r.pos == 0 ? 0 : r.pos - 1,
    static_cast<std::uint8_t>(r.mapq),
    r.cigar.size(),
    cigars.data(),
    r.rnext == "*" ? -1
                   : r.rnext == "=" ? tid
                                    : os.ref_table[r.rnext],
    r.pnext == 0 ? -1 : static_cast<std::int64_t>(r.pnext) - 1,
    r.tlen,
    r.seq.size(),
    r.seq.data(),
    quals,
    65535 // TODO: how to dynamicly adjust size with efficiency
  );*/
  auto aln = bam_init1();
  OBamStream::bam_set1(aln, r, tid, cigars, quals, mtid);
  if (auto c = sam_write1(os.bam, os.bam_header, aln); c < 0) { /* TODO: error handler */
    SPDLOG_ERROR("An error occured when writing to a BAM file");
  };
  return os;
}

}  // namespace biovoltron