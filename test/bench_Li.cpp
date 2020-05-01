#include "algorithm_327.h"
#include "algorithm_490.h"
#include "bench.hpp"
#include "cephes.h"
#include "hollik.h"
#include "Li2.hpp"
#include "Li3.hpp"
#include "Li4.hpp"
#include "Li5.hpp"
#include "Li6.hpp"
#include "Li.hpp"
#include "tsil_cpp.h"
#include <iostream>
#include <iomanip>
#include <gsl/gsl_sf_dilog.h>

extern "C" {
TSIL_REAL TSIL_dilog_real(TSIL_REAL x);
}

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

std::complex<double> hollik_Li2(std::complex<double> z) {
   double li2_re{}, li2_im{};
   hollik_dilog(std::real(z), std::imag(z), &li2_re, &li2_im);
   return { li2_re, li2_im };
}

long double tsil_Li2(long double x)
{
   return TSIL_dilog_real(x);
}

std::complex<long double> tsil_Li2(std::complex<long double> z)
{
   return TSIL_Dilog_(z);
}

std::complex<long double> tsil_Li3(std::complex<long double> z)
{
   return TSIL_Trilog_(z);
}

} // anonymous namespace

template <typename T, typename Fn>
void bench_fn(Fn f, const std::vector<T>& values, const std::string& name,
              const std::string& type)
{
   const auto total_time = polylogarithm::bench::time_in_seconds([&] {
         for (const auto& v: values) {
            polylogarithm::bench::do_not_optimize(f(v));
         }
      });

   std::cout << std::setw(24) << std::left << name << "type: " << std::setw(16)
             << std::left << type << "time: " << total_time << "s\n";
}

void print_line()
{
   std::cout << "----------------------------------------------------------------\n";
}

void print_headline(const std::string& text)
{
   print_line();
   std::cout << text << '\n';
   print_line();
}

int main() {
   using namespace polylogarithm::bench;

   const std::size_t N = 1000000;
   const auto min = -5.0;
   const auto max = 5.0;

   const auto values_d  = generate_random_scalars<double>(N, min, max);
   const auto values_l  = generate_random_scalars<long double>(N, min, max);
   const auto values_cd = generate_random_complexes<double>(N, min, max);
   const auto values_cl = generate_random_complexes<long double>(N, min, max);

   print_headline("Li2 (real)");

   bench_fn([&](double x) { return polylogarithm::Li2(x); }, values_d,
            "polylogarithm", "double");

   bench_fn([&](double x) { return gsl_Li2(x); }, values_d,
            "GSL", "double");

   bench_fn([&](double x) { return cephes_dilog(x); }, values_d,
            "cephes", "double");

   bench_fn([&](double x) { return algorithm_327(x); }, values_d,
            "algorithm 327", "double");

   bench_fn([&](double x) { return algorithm_490(x); }, values_d,
            "algorithm 490", "double");

   bench_fn([&](long double x) { return polylogarithm::Li2(x); }, values_l,
            "polylogarithm", "long double");

   bench_fn([&](long double x) { return tsil_Li2(x); }, values_l,
            "TSIL", "long double");

   print_headline("Li2 (complex)");

   bench_fn([&](std::complex<double> z) { return polylogarithm::Li2(z); },
            values_cd, "polylogarithm", "double");

   bench_fn([&](std::complex<double> z) { return gsl_Li2(z); },
            values_cd, "GSL", "double");

   bench_fn([&](std::complex<double> z) { return hollik_Li2(z); },
            values_cd, "Hollik", "double");

   bench_fn([&](std::complex<long double> z) { return polylogarithm::Li2(z); },
            values_cl, "polylogarithm",  "long double");

   bench_fn([&](std::complex<long double> z) { return tsil_Li2(z); },
            values_cl, "TSIL", "long double");

   print_headline("Li3 (complex)");

   bench_fn([&](std::complex<double> z) { return polylogarithm::Li3(z); },
            values_cd, "polylogarithm", "double");

   bench_fn([&](std::complex<long double> z) { return polylogarithm::Li3(z); },
            values_cl, "polylogarithm", "long double");

   bench_fn([&](std::complex<long double> z) { return tsil_Li3(z); },
            values_cl, "TSIL", "long double");

   print_headline("Li4 (complex)");

   bench_fn([&](std::complex<double> z) { return polylogarithm::Li4(z); },
            values_cd, "polylogarithm", "double");

   bench_fn([&](std::complex<long double> z) { return polylogarithm::Li4(z); },
            values_cd, "polylogarithm", "long double");

   print_headline("Li5 (complex)");

   bench_fn([&](std::complex<double> z) { return polylogarithm::Li5(z); },
            values_cd, "polylogarithm", "double");

   bench_fn([&](std::complex<long double> z) { return polylogarithm::Li5(z); },
            values_cd, "polylogarithm", "long double");

   print_headline("Li6 (complex)");

   bench_fn([&](std::complex<double> z) { return polylogarithm::Li6(z); },
            values_cd, "polylogarithm", "double");

   bench_fn([&](std::complex<long double> z) { return polylogarithm::Li6(z); },
            values_cd, "polylogarithm", "long double");

   // print_headline("Li(n,z) (complex)");

   // const auto values_cd_small = generate_random_complexes<double>(N/10, min, max);

   // bench_fn([&](std::complex<double> z) { return polylogarithm::Li(-100,z); },
   //          values_cd_small, "n=-100 polylogarithm", "double");

   // bench_fn([&](std::complex<double> z) { return polylogarithm::Li(-10,z); },
   //          values_cd_small, "n=-10  polylogarithm", "double");

   // bench_fn([&](std::complex<double> z) { return polylogarithm::Li(-6,z); },
   //          values_cd_small, "n=-6   polylogarithm", "double");

   // bench_fn([&](std::complex<double> z) { return polylogarithm::Li(6,z); },
   //          values_cd_small, "n=6    polylogarithm", "double");

   // bench_fn([&](std::complex<double> z) { return polylogarithm::Li(10,z); },
   //          values_cd_small, "n=10   polylogarithm", "double");

   // bench_fn([&](std::complex<double> z) { return polylogarithm::Li(100,z); },
   //          values_cd_small, "n=100  polylogarithm", "double");

   return 0;
}
