#include "bench.hpp"
#include "Li2.hpp"
#include "Li3.hpp"
#include "Li4.hpp"
#include "Li5.hpp"
#include "Li6.hpp"
#include "Li.hpp"
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

   std::cout << "Evaluation of real Li2 " << N << " times took: "
             << total_time << "s (GSL)\n";
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

   std::cout << "Evaluation of cmpl Li2 " << N << " times took: "
             << total_time << "s (GSL)\n";
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


void bench_Li3(std::size_t N, double min, double max)
{
   using namespace polylogarithm::bench;

   const auto values = generate_random_complexes(N, min, max);

   const double total_time = time_in_seconds([&values] {
         for (const auto& v: values) {
            volatile auto li = polylogarithm::Li3(v);
         }
      });

   std::cout << "Evaluation of cmpl Li3 " << N << " times took: "
             << total_time << "s\n";
}


void bench_Li4(std::size_t N, double min, double max)
{
   using namespace polylogarithm::bench;

   const auto values = generate_random_complexes(N, min, max);

   const double total_time = time_in_seconds([&values] {
         for (const auto& v: values) {
            volatile auto li = polylogarithm::Li4(v);
         }
      });

   std::cout << "Evaluation of cmpl Li4 " << N << " times took: "
             << total_time << "s\n";
}


void bench_Li5(std::size_t N, double min, double max)
{
   using namespace polylogarithm::bench;

   const auto values = generate_random_complexes(N, min, max);

   const double total_time = time_in_seconds([&values] {
         for (const auto& v: values) {
            volatile auto li = polylogarithm::Li5(v);
         }
      });

   std::cout << "Evaluation of cmpl Li5 " << N << " times took: "
             << total_time << "s\n";
}


void bench_Li6(std::size_t N, double min, double max)
{
   using namespace polylogarithm::bench;

   const auto values = generate_random_complexes(N, min, max);

   const double total_time = time_in_seconds([&values] {
         for (const auto& v: values) {
            volatile auto li = polylogarithm::Li6(v);
         }
      });

   std::cout << "Evaluation of cmpl Li6 " << N << " times took: "
             << total_time << "s\n";
}


void bench_Lin(long n, std::size_t N, double min, double max)
{
   using namespace polylogarithm::bench;

   const auto values = generate_random_complexes(N, min, max);

   const double total_time = time_in_seconds([&values,n] {
         for (const auto& v: values) {
            volatile auto li = polylogarithm::Li(n,v);
         }
      });

   std::cout << "Evaluation of cmpl Li(" << n << ") " << N << " times took: "
             << total_time << "s\n";
}


int main() {
   const std::size_t N = 1000000;
   const auto min = -10.;
   const auto max = 10.;

   bench_real_Li2(N, min, max);
   bench_cmpl_Li2(N, min, max);

   bench_real_Li2_GSL(N, min, max);
   bench_cmpl_Li2_GSL(N, min, max);

   bench_Li3(N, min, max);
   bench_Li4(N, min, max);
   bench_Li5(N, min, max);
   bench_Li6(N, min, max);

   bench_Lin(-100, N/10, min, max);
   bench_Lin(-10 , N/10, min, max);
   bench_Lin(-6  , N/10, min, max);
   bench_Lin( 6  , N/10, min, max);
   bench_Lin( 10 , N/10, min, max);
   bench_Lin( 100, N   , min, max);

   return 0;
}
