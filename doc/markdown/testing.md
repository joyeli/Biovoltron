# Testing Guide {#testing}

This guide outlines the standards and procedures for writing, running, and analyzing tests within the Biovoltron project. Adhering to these practices is crucial for maintaining code quality, stability, and maintainability.

[TOC]

## The Importance of Testing {#testing-importance}

Writing tests is a mandatory part of the development process. Every new feature, bug fix, or refactor must be accompanied by a comprehensive suite of tests that validate its correctness. Thoughtful, well-planned tests are required; trivial checks that do not meaningfully verify functionality are insufficient.

## Writing Effective Tests {#testing-writing-tests}

We use the [Catch2](https://github.com/catchorg/Catch2) framework for testing. To ensure clarity and organization, please follow these conventions.

### Test Structure

-   **TEST_CASE**: Each `TEST_CASE` should correspond to a single, logical unit of functionality (e.g., a specific function, a class method).
-   **SECTION**: Within a `TEST_CASE`, use `SECTION` blocks to test distinct sub-features, edge cases, or different input conditions. This isolates test logic and provides clearer output upon failure.

A good test validates behavior with meaningful assertions and comparisons against expected outcomes.

```cpp
#include "catch.hpp"
#include "biovoltron/my_feature.hpp"

TEST_CASE("MyFeature::doFeature - Validates core functionality", "[MyFeature]") {
  // Setup common to all sections can go here.

  SECTION("Test for positive input") {
    REQUIRE(my_feature::calculate(10) == 20);
  }

  SECTION("Test for negative input") {
    REQUIRE(my_feature::calculate(-5) == -10);
  }

  SECTION("Test for zero input") {
    REQUIRE(my_feature::calculate(0) == 0);
  }
}
```

### Test Case Naming Convention

To improve clarity and simplify test status analysis, all test cases must follow a standardized naming convention.

The format is `ClassName::MethodName - Behavior`, followed by a `[ClassName]` tag.

-   **Format**: `ClassName::MethodName - Behavior`
-   **Tag**: `[ClassName]`

The `ClassName` in both the name and the tag must refer to the primary class being tested, not any helper or secondary classes.

**Example:**

```cpp
TEST_CASE("MyCoolClass::doSomething - Handles positive numbers", "[MyCoolClass]") {
    // ... test logic ...
}
```

## Compiling for Coverage Analysis {#testing-compiling}

To generate code coverage data, you must enable the `GCOV` option during CMake configuration. This instruments the build to produce the necessary output files for coverage analysis.

```bash
# From your build directory
cmake -DGCOV=ON ..
make
```

## Running Tests {#testing-running-tests}

While the CI pipeline automatically runs tests, we strongly recommend running them locally before pushing changes to ensure they pass and to iterate more quickly.

### Local Test Execution

The main test executable is `biovoltron-test`, located in your build directory\'s `tests/` subdirectory.

1.  **Run all tests:**

```bash
# From your build directory
./tests/biovoltron-test
```

2.  **Run a specific Test Case:**

To save time, you can execute a single `TEST_CASE` by passing its name as an argument. This is highly efficient for focused development.

```bash
# From your build directory
./tests/biovoltron-test "MyFeature*"
```

## Generating Coverage Reports {#testing-coverage-reports}

After running the tests with a GCOV-enabled build, you can generate a detailed HTML coverage report.

### Requirements

-   **gcovr**: Version 8.0 or higher is required.

### Command

From your build directory, run the following command. The `--filter` option is crucial for targeting the report to the specific file(s) you have modified.

```bash
# From your build directory
gcovr \
  --html-details \
  --html-single-page \
  --gcov-ignore-parse-errors=suspicious_hits.warn_once_per_file \
  --exclude-unreachable-branches \
  --print-summary \
  -o coverage.html \
  --root ../include/biovoltron \
  --filter "../include/biovoltron/path/to/your/file.hpp" \
  .
```

This will generate a `coverage.html` file in your build directory. Open this file in a web browser to view the detailed analysis.

### Coverage Goals

The primary goal for all new contributions is to achieve **100% line-rate**.

While a 100% function-rate and branch-rate are also ideal and highly encouraged, the line-rate is the most critical metric we track.

## CI/CD Pipeline Integration {#testing-ci-cd}

The GitLab CI pipeline is configured to run all tests and generate a full coverage report. However, after pushing your changes, you must **manually initiate the pipeline runner** from the GitLab web interface to trigger the test suite.

The results will be available in the CI/CD job logs. This automated check serves as a final verification, but it should not replace local testing.
