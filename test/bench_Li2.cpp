#include "bench.hpp"
#include "Li2.hpp"
#include <iostream>


void bench_real_Li2(std::size_t N, double min, double max)
{
   using namespace polylogarithm::bench;

   const auto values = generate_random_doubles(N, min, max);

   const double total_time = time_in_milliseconds([&values] {
         for (const auto& v: values) {
            volatile auto li = polylogarithm::Li2(v);
         }
      });

   std::cout << "Evaluation of real Li2 " << N << " times took: "
             << total_time << "ms\n";
}


void bench_cmpl_Li2(std::size_t N, double min, double max)
{
   using namespace polylogarithm::bench;

   const auto values = generate_random_complexes(N, min, max);

   const double total_time = time_in_milliseconds([&values] {
         for (const auto& v: values) {
            volatile auto li = polylogarithm::Li2(v);
         }
      });

   std::cout << "Evaluation of cmpl Li2 " << N << " times took: "
             << total_time << "ms\n";
}


int main() {
   bench_real_Li2(1000000, -10., 10.);
   bench_cmpl_Li2(1000000, -10., 10.);

   return 0;
}
