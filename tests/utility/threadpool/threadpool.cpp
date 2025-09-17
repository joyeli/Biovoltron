#include <biovoltron/utility/threadpool/threadpool.hpp>
#include <chrono>
#include <thread>
#include <set>
#include <numeric>
#include <thread>
#include <catch.hpp>

TEST_CASE("ThreadPool::Execution - Executes tasks correctly", "[ThreadPool]") {

  SECTION("size check") {
    auto tp = biovoltron::make_threadpool();
    REQUIRE(tp.size() == std::thread::hardware_concurrency());
  }

  SECTION("sequential executation") {
    auto tp = biovoltron::make_threadpool(1);
    auto task = [](int x) {
      return x;
    };
    int n = 100;
    std::vector<int> v(n);
    std::vector<std::future<int>> res;

    iota(v.begin(), v.end(), 0);
    for (int i = 0; i < n; i++) {
      auto [id, fut] = tp.submit(task, i);
      res.emplace_back(std::move(fut));
    }

    std::vector<int> t;
    for (auto& f : res) {
      t.emplace_back(f.get());
    }
    REQUIRE(v == t);
  }

  SECTION("parallel executation") {
    int n = 1000;
    auto task = [](int x) {
      return 1LL * x * x;
    };

    std::vector<std::future<long long>> res;
    auto tp = biovoltron::make_threadpool();

    for (int i = 0; i < n; i++) {
      res.emplace_back(tp.submit(task, i).second);
    }

    std::set<long long> st;
    for (auto& f : res) {
      st.insert(f.get());
    }
    REQUIRE(st.size() == n);
  }
}



