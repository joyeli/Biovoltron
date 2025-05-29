#pragma once

#include <ostream>

namespace biovoltron {
using Genotype = std::pair<std::uint8_t, std::uint8_t>;
}

namespace std {
inline auto&
operator<<(ostream& os, biovoltron::Genotype genotype) {
  return os << +genotype.first << "|" << +genotype.second;
}
}  // namespace std