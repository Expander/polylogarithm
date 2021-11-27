#include "alt.h"
#include "bench.hpp"
#include "c_wrappers.h"
#include "fortran_wrappers.h"
#include "Cl2.hpp"
#include "Cl3.hpp"
#include "Cl4.hpp"
#include "Cl5.hpp"
#include "Cl6.hpp"
#include "Li2.hpp"
#include "Li3.hpp"
#include "Li4.hpp"
#include <iostream>
#include <iomanip>

#ifdef ENABLE_GSL
#include <gsl/gsl_sf_clausen.h>
#endif

namespace {

#ifdef ENABLE_FORTRAN

double poly_Cl2_fortran(double x) {
   double res{};
   cl2_fortran(&x, &res);
   return res;
}

#endif

double Cl2_via_Li2(double x) noexcept
{
   return std::imag(polylogarithm::Li2(std::polar(1.0, x)));
}

double Cl3_via_Li3(double x) noexcept
{
   return std::real(polylogarithm::Li3(std::polar(1.0, x)));
}

double Cl4_via_Li4(double x) noexcept
{
   return std::imag(polylogarithm::Li4(std::polar(1.0, x)));
}

long double Cl2_via_Li2(long double x) noexcept
{
   return std::imag(polylogarithm::Li2(std::polar(1.0L, x)));
}

long double Cl3_via_Li3(long double x) noexcept
{
   return std::real(polylogarithm::Li3(std::polar(1.0L, x)));
}

long double Cl4_via_Li4(long double x) noexcept
{
   return std::imag(polylogarithm::Li4(std::polar(1.0L, x)));
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

void print_headline_2(const std::string& text)
{
   print_line('-');
   std::cout << text << '\n';
   print_line('-');
}

template<typename T, typename U>
void bench(const T& values_d, const U& values_l)
{
   print_headline_2("Cl2");

   bench_fn([&](double x) { return polylogarithm::Cl2(x); }, values_d,
            "polylogarithm C++", "double");

   bench_fn([&](double x) { return cl2(x); }, values_d,
            "polylogarithm C", "double");

#ifdef ENABLE_FORTRAN
   bench_fn([&](double x) { return poly_Cl2_fortran(x); }, values_d,
            "polylogarithm Fortran", "double");
#endif

   bench_fn([&](double x) { return Cl2_via_Li2(x); }, values_d,
            "via Li2 C++", "double");

   bench_fn([&](double x) { return clausen_2_bernoulli(x); }, values_d,
            "Bernoulli", "double");

#ifdef ENABLE_GSL
   bench_fn([&](double x) { return gsl_sf_clausen(x); }, values_d,
            "GSL", "double");
#endif

   bench_fn([&](double x) { return clausen_2_koelbig(x); }, values_d,
            "Koelbig C", "double");

   bench_fn([&](double x) { return clausen_2_pade(x); }, values_d,
            "Pade C", "double");

   bench_fn([&](double x) { return clausen_2_wu(x); }, values_d,
            "Wu", "double");

   bench_fn([&](long double x) { return polylogarithm::Cl2(x); }, values_l,
            "polylogarithm C++", "long double");

   bench_fn([&](long double x) { return cl2l(x); }, values_l,
            "polylogarithm C", "long double");

   bench_fn([&](long double x) { return clausen_2l_koelbig(x); }, values_l,
            "Koelbig C", "long double");

   bench_fn([&](long double x) { return Cl2_via_Li2(x); }, values_l,
            "via Li2 C++", "long double");

   print_headline_2("Cl3");

   bench_fn([&](double x) { return polylogarithm::Cl3(x); }, values_d,
            "polylogarithm C++", "double");

   bench_fn([&](double x) { return cl3(x); }, values_d,
            "polylogarithm C", "double");

   bench_fn([&](double x) { return Cl3_via_Li3(x); }, values_d,
            "via Li3 C++", "double");

   bench_fn([&](double x) { return clausen_3_pade(x); }, values_d,
            "Pade", "double");

   bench_fn([&](double x) { return clausen_3_wu(x); }, values_d,
            "Wu", "double");

   bench_fn([&](long double x) { return polylogarithm::Cl3(x); }, values_l,
            "polylogarithm C++", "long double");

   bench_fn([&](long double x) { return Cl3_via_Li3(x); }, values_l,
            "via Li3 C++", "long double");

   print_headline_2("Cl4");

   bench_fn([&](double x) { return polylogarithm::Cl4(x); }, values_d,
            "polylogarithm C++", "double");

   bench_fn([&](double x) { return cl4(x); }, values_d,
            "polylogarithm C", "double");

   bench_fn([&](double x) { return Cl4_via_Li4(x); }, values_d,
            "via Li4 C++", "double");

   bench_fn([&](long double x) { return polylogarithm::Cl4(x); }, values_l,
            "polylogarithm C++", "long double");

   bench_fn([&](long double x) { return Cl4_via_Li4(x); }, values_l,
            "via Li4 C++", "long double");

   print_headline_2("Cl5");

   bench_fn([&](double x) { return polylogarithm::Cl5(x); }, values_d,
            "polylogarithm C++", "double");

   bench_fn([&](long double x) { return polylogarithm::Cl5(x); }, values_l,
            "polylogarithm C++", "long double");

   print_headline_2("Cl6");

   bench_fn([&](double x) { return polylogarithm::Cl6(x); }, values_d,
            "polylogarithm C++", "double");

   bench_fn([&](long double x) { return polylogarithm::Cl6(x); }, values_l,
            "polylogarithm C++", "long double");
}

int main()
{
   using polylogarithm::bench::generate_random_scalars;
   using polylogarithm::bench::generate_random_complexes;

   const std::size_t N = 1000000;
   const auto pi = 3.1415926535897932;
   const auto min = -8*pi;
   const auto max = 8*pi;

   // range [0,pi), where no range reduction is necessary
   const auto values_d_redu = generate_random_scalars<double>(N, 0, pi);
   const auto values_l_redu = generate_random_scalars<long double>(N, 0, pi);
   // extended range, where range reduction is necessary
   const auto values_d_full = generate_random_scalars<double>(N, min, max);
   const auto values_l_full = generate_random_scalars<long double>(N, min, max);

   print_headline_1("Benchmark without range reduction");
   bench(values_d_redu, values_l_redu);

   print_headline_1("Benchmark with range reduction");
   bench(values_d_full, values_l_full);

   return 0;
}
