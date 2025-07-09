#pragma once

#include <biovoltron/utility/archive/gzstream.hpp>
#include <spdlog/spdlog.h>
#include <filesystem>

namespace fs = std::filesystem;

namespace biovoltron::detail {

enum FILE_FORMAT {
  FASTA = 1 << 0,
  FASTQ = 1 << 1,
  SAM   = 1 << 2,
  BAM   = 1 << 3,
  GZ    = 1 << 4,
  ERROR = 1 << 5
};

auto parse_file_format(fs::path p) {
  auto ext = p.extension();
  int type = 0;
  if (ext == ".gz") {
    type |= FILE_FORMAT::GZ;
    p.replace_extension("");
  }
  if (ext == ".fa" || ext == ".fasta") {
    type |= FILE_FORMAT::FASTA;
  } else if (ext == ".fq" || ext == ".fastq") {
    type |= FILE_FORMAT::FASTQ;
  } else if (ext == ".sam") {
    type |= FILE_FORMAT::SAM;
  } else if (ext == ".bam") {
    type |= FILE_FORMAT::BAM;
  } else {
    type |= FILE_FORMAT::ERROR;
  }
  if (type & (FILE_FORMAT::GZ | FILE_FORMAT::BAM)) {
    type |= FILE_FORMAT::ERROR;
  }
  return type;
}


auto check_path_exists(const std::filesystem::path& path) {
  if (!std::filesystem::exists(path)) {
    SPDLOG_ERROR("{} doesn't exists", path.string());
    throw std::invalid_argument("File not found");
  }
}

template<class IStream = std::ifstream> 
  requires std::is_base_of_v<std::ifstream, IStream>
auto open_input_file(const std::filesystem::path& path) {
  check_path_exists(path);
  IStream fin(std::filesystem::absolute(path));
  if (!fin.is_open()) {
    SPDLOG_ERROR("Cannot open {}", path.string());
    throw std::invalid_argument("Cannot open file");
  }
  return fin;
};

template<class OStream = std::ofstream> 
  requires std::is_base_of_v<std::ofstream, OStream>
auto open_output_file(const std::filesystem::path& path) {
  OStream fout(std::filesystem::absolute(path));
  if (!fout.is_open()) {
    SPDLOG_ERROR("Cannot open {}", path.string());
    throw std::invalid_argument("Cannot open file");
  }
  return fout;
}

} // namespace biovoltron::detail