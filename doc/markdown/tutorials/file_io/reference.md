# Working with ReferenceRecord {#tutorials-fileio-reference}

[TOC]

In this tutorial, we will learn how to work with `biovoltron::ReferenceRecord`, a data structure designed to store reference genomes in Biovoltron.
A `ReferenceRecord` can be constructed directly from FASTA files, saved into a compact binary format, and later reloaded for efficient reuse.

## Step 1. Include the Necessary Headers

To use `ReferenceRecord`, you need to include the following header:

```cpp
#include <biovoltron/file_io/reference.hpp>
```

## Step 2. Creating a ReferenceRecord Object

You can create a `ReferenceRecord` without encoding:

```cpp
// Sequence stored as plain std::string
biovoltron::ReferenceRecord<false> ref;
ref.seq = "ACTGAACAGTCCATCGTGACTGGACTGACTGA";
```

Or with encoding:

```cpp
biovoltron::ReferenceRecord<true> ref_encoded;
ref_encoded.seq = biovoltron::istring{0, 1, 2, 3};
```

* `ReferenceRecord<false>` keeps bases as normal characters (`A`, `C`, `G`, `T`, `N`).
* `ReferenceRecord<true>` stores bases as integers (0–3).

## Step 3. Loading from a FASTA File

`ReferenceRecord` integrates seamlessly with FASTA input.
You can parse a FASTA file directly into a `ReferenceRecord`.

> Note: The following example assumes you have a FASTA file named `GRCh38` in your working directory.

```cpp
#include <fstream>
#include <iostream>
#include <biovoltron/file_io/reference.hpp>

using namespace biovoltron;

int main() {
  // Open a FASTA reference genome file
  auto fin = std::ifstream{"GRCh38"};

  // Create an empty ReferenceRecord
  auto ref = ReferenceRecord<false>{ .species = "Human" };

  // Parse the FASTA data into ref
  fin >> ref;

  // Print some basic information
  std::cout << "Species: " << ref.species << "\n";
  std::cout << "Chromosomes: " << ref.chr_num << "\n";

  return 0;
}
```


This automatically:

* Counts bases (`A`, `C`, `G`, `T`, `N`).
* Records chromosome names and cumulative end positions.
* Detects unknown (`N`) regions and stores them in `unknown_intervals`.

## Step 4. Inspecting the Reference

Once loaded, you can inspect different parts of the reference:

```cpp
// Chromosome names
for (const auto& name : ref.chr_names) {
  std::cout << "Chromosome: " << name << "\n";
}

// Base counts (0:A, 1:C, 2:G, 3:T, 4:N)
for (size_t i = 0; i < ref.base_cnt.size(); ++i) {
  std::cout << "Base[" << i << "] = " << ref.base_cnt[i] << "\n";
}

// End positions of chromosomes
for (const auto& pos : ref.chr_end_pos) {
  std::cout << "Chr end at: " << pos << "\n";
}

// Unknown intervals
for (const auto& interval : ref.unknown_intervals) {
  std::cout << "Unknown region: " 
            << interval[0] << " - " << interval[1] << "\n";
}
```

## Step 5. Recovering Original Sequence

If you want to reconstruct the original sequence **with `N` in unknown regions**, use `origin_seq()`:

```cpp
auto original_sequence = ref.origin_seq();
std::cout << "Recovered sequence = " << original_sequence << "\n";
```

This restores the `N` bases in all unknown intervals.

## Step 6. Saving and Loading in Binary Format

`ReferenceRecord` can be saved into a compact `.bfa` binary file and later loaded back.

> Note: This example assumes you have a FASTA file named `GRCh38` to create the initial record.

```cpp
#include <fstream>
#include <iostream>
#include <cassert>
#include <biovoltron/file_io/reference.hpp>

using namespace biovoltron;

int main() {
  // First, create and load a ReferenceRecord from a FASTA file.
  auto ref = ReferenceRecord<false>{};
  {
    auto fin = std::ifstream{"GRCh38"};
    if (!fin) {
      std::cerr << "Error: Cannot open GRCh38. Please create this file.\n";
      // As a fallback for the example, create a dummy record.
      ref.seq = "ACGTN";
      ref.chr_names.push_back("chr1");
      ref.chr_end_pos.push_back(5);
      ref.base_cnt = {1,1,1,1,1};
      ref.chr_num = 1;
      ref.species = "Dummy";
    } else {
      fin >> ref;
    }
  }

  // Save the record to a binary .bfa file
  {
    auto fout = std::ofstream{"GRCh38.bfa", std::ios::binary};
    ref.save(fout);
    std::cout << "Saved reference to GRCh38.bfa" << std::endl;
  }

  // Load the record from the binary .bfa file
  auto loaded_ref = ReferenceRecord<false>{};
  {
    auto fin = std::ifstream{"GRCh38.bfa", std::ios::binary};
    loaded_ref.load(fin);
    std::cout << "Loaded reference from GRCh38.bfa" << std::endl;
  }

  // Verify that the saved and loaded records are identical
  assert(ref == loaded_ref);
  std::cout << "Verification successful: original and loaded records are identical." << std::endl;

  return 0;
}
```


<span class="next_section_button">
[Read Next: ReferenceRecord Modules](structbiovoltron_1_1ReferenceRecord.html)
</span>
