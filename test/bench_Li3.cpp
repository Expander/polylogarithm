#include "bench.hpp"
#include "Li3.hpp"
#include <iostream>


void bench_Li3(std::size_t N, double min, double max)
{
   using namespace polylogarithm::bench;

   const auto values = generate_random_complexes(N, min, max);

   const double total_time = time_in_seconds([&values] {
         for (const auto& v: values) {
            volatile auto li = polylogarithm::Li3(v);
         }
      });

   std::cout << "Evaluation of complex Li3 " << N << " times took: "
             << total_time << "s\n";
}


int main() {
   bench_Li3(1000000, -10., 10.);

   return 0;
}
