# README: WFPOA Compilation and Execution Guide

WFPOA (Wavefront Partial Order Alignment) is a tool designed for sequence alignment and consensus generation. This guide explains how to compile and execute WFPOA using different methods.

## Sample Usage
```cpp
auto alignment_engine = biovoltron::SimdAlignmentEngine::Create(
    biovoltron::AlignmentType::kNW, 0, -1, -1);
biovoltron::Graph graph{};

// sequence to graph alignment
for (const auto& it : sequences) {
    auto alignment = alignment_engine->Align(it, graph, 0);
    graph.AddAlignment(alignment, it);
}

auto c = graph.GenerateConsensus();
```

## Compilation
Compile WFPOA using one of the following methods:

### Default Method (`wf_blunt`)
```sh
g++ test.cpp -std=c++17 -mavx2 -O3 -o wfpoaBlunt 
```

### `wf_unit` Method
```sh
g++ -D WFUNIT test.cpp -std=c++17 -mavx2 -O3 -o wfpoaUnit 
```

### `wf_arrow` Method
```sh
g++ test.cpp -DWFARROW -std=c++17 -mavx2 -O3 -o wfpoaArrow 
```

## Execution
Run WFPOA using the following command:
```sh
./wfpoa{method} {cut_threshold} {filepath}
```

### Parameters:
- `{method}`: The method used (`Blunt`, `Unit`, or `Arrow`).
- `{cut_threshold}`: A numerical threshold for processing (default: `10`, conservative setting: `20`).
- `{filepath}`: The path to the input sequence file.

### Example Execution
```sh
./wfpoaBlunt 10 f1.fa
```

For CCS reads, a `cut_threshold` of `10` is usually sufficient, but setting it to `20` ensures a more conservative approach.

---

This README provides a concise explanation of WFPOA's usage. Adjust parameters based on your specific needs.
