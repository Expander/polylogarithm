#include "bench.hpp"
#include "Li2.hpp"
#include <iostream>


void bench_real_Li2(std::size_t N, double min, double max)
{
   using namespace polylogarithm::bench;

   const auto values = generate_random_doubles(N, min, max);
   double total_time = 0.;

   for (const auto& v: values) {
      total_time += time([v] { return polylogarithm::Li2(v); });
   }

   std::cout << "Evaluation of real Li2 " << N << " times took: "
             << total_time*1000. << "ms\n";
}


void bench_cmpl_Li2(std::size_t N, double min, double max)
{
   using namespace polylogarithm::bench;

   const auto values = generate_random_complexes(N, min, max);
   double total_time = 0.;

   for (const auto& v: values) {
      total_time += time([v] { return polylogarithm::Li2(v); });
   }

   std::cout << "Evaluation of complex Li2 " << N << " times took: "
             << total_time*1000. << "ms\n";
}


int main() {
   bench_real_Li2(1000000, -10., 10.);
   bench_cmpl_Li2(1000000, -10., 10.);

   return 0;
}
