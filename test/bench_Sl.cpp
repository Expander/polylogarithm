#include "bench.hpp"
#include "Sl.hpp"
#include <iostream>
#include <iomanip>
#include <string>

template <typename T, typename Fn>
void bench_fn(Fn f, const std::vector<T>& values, const std::string& name,
              const std::string& type)
{
   // warm-up
   for (const auto& v: values) {
      polylogarithm::bench::do_not_optimize(f(v));
   }

   const auto total_time = polylogarithm::bench::time_in_seconds([&] {
         for (const auto& v: values) {
            polylogarithm::bench::do_not_optimize(f(v));
         }
      });

   std::cout << std::setw(24) << std::left << name << "type: " << std::setw(16)
             << std::left << type << "time: " << total_time << "s\n";
}

void print_line(char c)
{
   for (int i = 0; i < 60; ++i) {
      std::cout << c;
   }
   std::cout << '\n';
}

void print_headline_1(const std::string& text)
{
   print_line('=');
   std::cout << text << '\n';
   print_line('=');
}

template<typename T>
void bench(const T& values_d)
{
   const int ni[] = {
       1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
      11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
      21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
      31, 1000, 10000, 100000, 1000000, 10000000
   };

   for (const auto n: ni) {
      bench_fn([&](double x) { return polylogarithm::Sl(n,x); }, values_d,
               std::string("Sl(") + std::to_string(n) + ", x)", "double");
   }
}

int main()
{
   using polylogarithm::bench::generate_random_scalars;

   const std::size_t N = 10000000;
   const auto pi = 3.1415926535897932;
   const auto min = -8*pi;
   const auto max = 8*pi;

   // range [0,pi), where no range reduction is necessary
   const auto values_d_redu = generate_random_scalars<double>(N, 0, pi);
   // extended range, where range reduction is necessary
   const auto values_d_full = generate_random_scalars<double>(N, min, max);

   print_headline_1("Benchmark without range reduction");
   bench(values_d_redu);

   print_headline_1("Benchmark with range reduction");
   bench(values_d_full);

   return 0;
}
