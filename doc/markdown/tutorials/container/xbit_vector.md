# Working with XbitVector {#tutorials-container-xbitvector}

\[TOC]

In this tutorial, we will explore how to effectively use **`biovoltron::DibitVector`** and **`biovoltron::QuadbitVector`**, two space-efficient containers from Biovoltron.
They are based on a general template `XbitVector<N, Block, Allocator>` and work just like `std::vector`, but pack multiple values into fewer bits.

* **`DibitVector`** packs **2 bits per value** (4 possible values: 0, 1, 2, 3).
  → Useful for DNA bases (`A, C, G, T`).
* **`QuadbitVector`** packs **4 bits per value** (16 possible values: 0–15).

## Step 1. Include the Necessary Header

To use these containers, include the following header:

```cpp
#include <biovoltron/container/xbit_vector.hpp>
```

## Step 2. Creating a DibitVector

You can create a `DibitVector` just like you would create a `std::vector`:

```cpp
#include <cassert>
#include <iostream>
#include <biovoltron/container/xbit_vector.hpp>

int main() {
  // Create a DibitVector with 5 elements, all initialized to 2
  biovoltron::DibitVector<> v(5, 2);

  std::cout << "First element: " << +v.front() << "\n";
}
```

Here, `+v.front()` prints `2`.
Notice the `+` before `v.front()` — since the value type is `uint8_t`,
we cast to `int` for printing.

## Step 3. Using Initializer Lists

You can directly initialize with values:

```cpp
biovoltron::DibitVector<> v = {0, 1, 2, 3, 2, 1};
```

This creates a vector of six elements: `[0,1,2,3,2,1]`.

## Step 4. Iterating

```cpp
#include <iostream>
#include <biovoltron/container/xbit_vector.hpp>
#include <algorithm>

int main() {
  biovoltron::DibitVector<> v = {3, 2, 1, 0};

  // Range-based for loop
  for (auto x : v)
    std::cout << +x << " ";  // prints: 3 2 1 0
}
```

## Step 5. Specialized Operations

`XbitVector` has some **extra utilities**, such as **`flip()`** → flip all values at once

Example:

```cpp
biovoltron::DibitVector<> v = {0, 1, 2, 3};
v.flip();

for (auto x : v)
  std::cout << +x << " ";  // prints: 3 2 1 0
```

## Step 6. DNA Example with DibitVector

A common use case is encoding DNA bases:

```cpp
biovoltron::DibitVector<> dna = {0, 1, 2, 3}; // A=0, C=1, G=2, T=3

auto base_view = std::views::transform([](auto c){ return "ACGT"[c]; });
std::ranges::copy(dna | base_view,
                  std::ostream_iterator<char>{std::cout, ""});
// Output: ACGT
```

You can also compute the reverse complement:

```cpp
auto comp_view = std::views::transform([](auto c){ return 0b11u - c; });
for (auto c : dna | std::views::reverse | comp_view | base_view)
  std::cout << c;  // Output: ACGT -> reverse complement -> ACGT
```

## Step 7. QuadbitVector Example

`QuadbitVector` works the same, but allows values from 0 to 15:

```cpp
biovoltron::QuadbitVector<> q = {1, 5, 9, 15};
for (auto x : q)
  std::cout << +x << " ";  // prints: 1 5 9 15
```

<span class="next_section_button">
[Read Next: XbitVector Modules](classbiovoltron_1_1detail_1_1XbitVector.html)
</span>
