#pragma once

#include <biovoltron/file_io/core/header.hpp>
#include <biovoltron/file_io/core/record.hpp>

namespace biovoltron {

/**
 * @ingroup file_io
 * @brief TBA
 */
struct WigHeader : Header {
  constexpr static auto START_SYMBOLS =
    std::array{"browser", "#", "track", "variableStep", "fixedStep"};
};

/**
 * @ingroup file_io
 * @brief TBA
 */
struct WigVarStepRecord : HeaderableRecord {
  WigHeader* header = nullptr;

  std::uint32_t start{};
  float value{};
};

/**
 * @ingroup file_io
 * @brief TBA
 */
struct WigFixedStepRecord : HeaderableRecord {
  WigHeader* header = nullptr;

  float value{};
};

}  // namespace biovoltron
