#include "alt.h"
#include "bench.hpp"
#include "Cl2.hpp"
#include "Cl3.hpp"
#include "Cl4.hpp"
#include "Cl5.hpp"
#include "Cl6.hpp"
#include "Li2.hpp"
#include <iostream>
#include <iomanip>

#ifdef ENABLE_GSL
#include <gsl/gsl_sf_clausen.h>
#endif

namespace {

double Cl2_via_Li2(double x) noexcept
{
   return std::imag(polylogarithm::Li2(std::exp(std::complex<double>(0.0, x))));
}

long double Cl2_via_Li2(long double x) noexcept
{
   return std::imag(polylogarithm::Li2(std::exp(std::complex<long double>(0.0L, x))));
}

} // anonymous namespace

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
   using polylogarithm::bench::generate_random_scalars;
   using polylogarithm::bench::generate_random_complexes;

   const std::size_t N = 1000000;
   const auto pi = 3.1415926535897932;
   const auto min = -8*pi;
   const auto max = 8*pi;

   const auto values_d  = generate_random_scalars<double>(N, min, max);
   const auto values_l  = generate_random_scalars<long double>(N, min, max);

   print_headline("Cl2");

   bench_fn([&](double x) { return polylogarithm::Cl2(x); }, values_d,
            "polylogarithm C++", "double");

   bench_fn([&](double x) { return Cl2_via_Li2(x); }, values_d,
            "via Li2 C++", "double");

   bench_fn([&](double x) { return koelbig_cl2(x); }, values_d,
            "Koelbig C", "double");

#ifdef ENABLE_GSL
   bench_fn([&](double x) { return gsl_sf_clausen(x); }, values_d,
            "GSL", "double");
#endif

   bench_fn([&](double x) { return clausen_2_wu(x); }, values_d,
            "Wu", "double");

   bench_fn([&](long double x) { return polylogarithm::Cl2(x); }, values_l,
            "polylogarithm C++", "long double");

   bench_fn([&](long double x) { return Cl2_via_Li2(x); }, values_l,
            "via Li2 C++", "long double");

   print_headline("Cl3");

   bench_fn([&](double x) { return polylogarithm::Cl3(x); }, values_d,
            "polylogarithm C++", "double");

   bench_fn([&](double x) { return clausen_3_wu(x); }, values_d,
            "Wu", "double");

   bench_fn([&](long double x) { return polylogarithm::Cl3(x); }, values_l,
            "polylogarithm C++", "long double");

   print_headline("Cl4");

   bench_fn([&](double x) { return polylogarithm::Cl4(x); }, values_d,
            "polylogarithm C++", "double");

   bench_fn([&](long double x) { return polylogarithm::Cl4(x); }, values_l,
            "polylogarithm C++", "long double");

   print_headline("Cl5");

   bench_fn([&](double x) { return polylogarithm::Cl5(x); }, values_d,
            "polylogarithm C++", "double");

   bench_fn([&](long double x) { return polylogarithm::Cl5(x); }, values_l,
            "polylogarithm C++", "long double");

   print_headline("Cl6");

   bench_fn([&](double x) { return polylogarithm::Cl6(x); }, values_d,
            "polylogarithm C++", "double");

   bench_fn([&](long double x) { return polylogarithm::Cl6(x); }, values_l,
            "polylogarithm C++", "long double");

   return 0;
}
