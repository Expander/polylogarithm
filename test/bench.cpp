#include "Li2.hpp"
#include <benchmark/benchmark.h>

template <class ...ExtraArgs>
static void BM_Li2_real_poly_cpp(benchmark::State& state, ExtraArgs&&... extra_args)
{
  for (auto _ : state) {
     benchmark::DoNotOptimize(polylogarithm::Li2(extra_args...));
  }
}

BENCHMARK_CAPTURE(BM_Li2_real_poly_cpp, double_test_unity, 1.1);

BENCHMARK_MAIN();
