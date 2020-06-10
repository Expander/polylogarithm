#include "alt.h"
#include "Li2.hpp"
#include <benchmark/benchmark.h>

static const double points[] = {
   +2.1, +2.0, +1.5, +1.0, +0.6, +0.5, +0.1, +0.0,
   -2.1, -2.0, -1.5, -1.0, -0.6, -0.5, -0.1, -0.0
};

static void BM_Li2_poly_cpp(benchmark::State& state)
{
  for (auto _ : state) {
     benchmark::DoNotOptimize(polylogarithm::Li2(points[state.range(0)]));
  }
}

static void BM_Li2_poly_c(benchmark::State& state)
{
  for (auto _ : state) {
     benchmark::DoNotOptimize(li2(points[state.range(0)]));
  }
}

BENCHMARK(BM_Li2_poly_cpp)->DenseRange(0, (sizeof(points)/sizeof(points[0])), 1);
BENCHMARK(BM_Li2_poly_c)  ->DenseRange(0, (sizeof(points)/sizeof(points[0])), 1);

BENCHMARK_MAIN();
