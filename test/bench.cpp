#include "Li2.hpp"
#include <benchmark/benchmark.h>

template <class ...ExtraArgs>
static void BM_Li2_poly_cpp(benchmark::State& state, ExtraArgs&&... extra_args)
{
  for (auto _ : state) {
     benchmark::DoNotOptimize(polylogarithm::Li2(extra_args...));
  }
}


BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_unity,  0.0);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_unity,  0.1);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_unity,  0.4);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_unity,  0.5);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_unity,  0.6);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_unity,  0.9);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_unity,  1.0);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_unity,  1.2);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_unity,  1.5);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_unity,  1.9);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_unity,  2.0);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_unity,  2.1);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_unity,  3.0);

BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_unity, -0.0);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_unity, -0.1);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_unity, -0.4);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_unity, -0.5);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_unity, -0.6);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_unity, -0.9);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_unity, -1.0);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_unity, -1.2);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_unity, -1.5);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_unity, -1.9);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_unity, -2.0);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_unity, -2.1);
BENCHMARK_CAPTURE(BM_Li2_poly_cpp, double_unity, -3.0);

BENCHMARK_MAIN();
