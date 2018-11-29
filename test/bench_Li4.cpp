#include "bench.hpp"
#include "Li4.hpp"
#include <iostream>


void bench_Li4(std::size_t N, double min, double max)
{
   using namespace polylogarithm::bench;

   const auto values = generate_random_complexes(N, min, max);

   const double total_time = time_in_seconds([&values] {
         for (const auto& v: values) {
            volatile auto li = polylogarithm::Li4(v);
         }
      });

   std::cout << "Evaluation of complex Li4 " << N << " times took: "
             << total_time << "s\n";
}


int main() {
   bench_Li4(1000000, -10., 10.);

   return 0;
}
