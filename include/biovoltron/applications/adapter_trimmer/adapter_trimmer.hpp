#pragma once

#include <concepts>
#include <biovoltron/applications/adapter_trimmer/paired_end/trimmer.hpp>
#include <biovoltron/applications/adapter_trimmer/single_end/trimmer.hpp>

namespace biovoltron {

namespace earrings {
  // TODO: make a concept for acceptable record
}

template<class Record>
  requires std::derived_from<Record, FastaRecord<Record::encoded>>
  // TODO: || std::same_as<Record, BamRecord>
using PairedEndAdapterTrimmer = paired::AdapterTrimmer<Record>;

template<class Record>
  requires std::derived_from<Record, FastaRecord<Record::encoded>>
// TODO: || std::same_as<Record, BamRecord>
using SingleEndAdapterTrimmer = single::AdapterTrimmer<Record>;

} // namespace biovoltron
