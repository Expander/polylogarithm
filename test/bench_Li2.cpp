#include "bench.hpp"
#include "Li2.hpp"
#include <iostream>
#include <gsl/gsl_sf_dilog.h>

namespace {

double gsl_Li2(double x)
{
   gsl_sf_result li2_gsl{};
   gsl_sf_dilog_e(x, &li2_gsl);
   return li2_gsl.val;
}


std::complex<double> gsl_Li2(std::complex<double> z)
{
   gsl_sf_result li2_gsl_re{}, li2_gsl_im{};
   gsl_sf_complex_dilog_e(std::abs(z), std::arg(z), &li2_gsl_re, &li2_gsl_im);
   return {li2_gsl_re.val, li2_gsl_im.val};
}

} // anonymous namespace


void bench_real_Li2_GSL(std::size_t N, double min, double max)
{
   using namespace polylogarithm::bench;

   const auto values = generate_random_doubles(N, min, max);

   const double total_time = time_in_seconds([&values] {
         for (const auto& v: values) {
            volatile auto li = gsl_Li2(v);
         }
      });

   std::cout << "Evaluation of real Li2 (via GSL) " << N << " times took: "
             << total_time << "s\n";
}


void bench_real_Li2(std::size_t N, double min, double max)
{
   using namespace polylogarithm::bench;

   const auto values = generate_random_doubles(N, min, max);

   const double total_time = time_in_seconds([&values] {
         for (const auto& v: values) {
            volatile auto li = polylogarithm::Li2(v);
         }
      });

   std::cout << "Evaluation of real Li2 " << N << " times took: "
             << total_time << "s\n";
}


void bench_cmpl_Li2_GSL(std::size_t N, double min, double max)
{
   using namespace polylogarithm::bench;

   const auto values = generate_random_complexes(N, min, max);

   const double total_time = time_in_seconds([&values] {
         for (const auto& v: values) {
            volatile auto li = gsl_Li2(v);
         }
      });

   std::cout << "Evaluation of cmpl Li2 (via GSL) " << N << " times took: "
             << total_time << "s\n";
}


void bench_cmpl_Li2(std::size_t N, double min, double max)
{
   using namespace polylogarithm::bench;

   const auto values = generate_random_complexes(N, min, max);

   const double total_time = time_in_seconds([&values] {
         for (const auto& v: values) {
            volatile auto li = polylogarithm::Li2(v);
         }
      });

   std::cout << "Evaluation of cmpl Li2 " << N << " times took: "
             << total_time << "s\n";
}


int main() {
   bench_real_Li2(1000000, -10., 10.);
   bench_cmpl_Li2(1000000, -10., 10.);

   bench_real_Li2_GSL(1000000, -10., 10.);
   bench_cmpl_Li2_GSL(1000000, -10., 10.);

   return 0;
}
