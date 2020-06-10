#include "Li2.hpp"
#include <benchmark/benchmark.h>

static void BM_Li2_real_poly_cpp(benchmark::State& state)
{
  for (auto _ : state) {
     benchmark::DoNotOptimize(polylogarithm::Li2(1.1));
  }
}

BENCHMARK(BM_Li2_real_poly_cpp);

BENCHMARK_MAIN();
