#include "Li2.hpp"
#include <benchmark/benchmark.h>

template <class ...ExtraArgs>
static void BM_Li2_poly_cpp(benchmark::State& state, ExtraArgs&&... extra_args)
{
  for (auto _ : state) {
     benchmark::DoNotOptimize(polylogarithm::Li2(extra_args...));
  }
}

BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_p2_1,  2.1);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_p2_0,  2.0);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_p1_5,  1.5);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_p1_0,  1.0);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_p0_6,  0.6);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_p0_5,  0.5);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_p0_1,  0.1);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_p0_0,  0.0);

BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_m0_0, -0.0);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_m0_1, -0.1);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_m0_5, -0.5);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_m0_6, -0.6);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_m1_0, -1.0);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_m1_5, -1.5);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_m2_0, -2.0);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_m2_1, -2.1);

BENCHMARK_MAIN();
