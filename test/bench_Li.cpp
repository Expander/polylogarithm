#include "alt.h"
#include "bench.hpp"
#include "c_wrappers.h"
#include "fortran_wrappers.h"
#include "Li2.hpp"
#include "Li3.hpp"
#include "Li4.hpp"
#include "Li5.hpp"
#include "Li6.hpp"
#include "Li.hpp"
#include <iostream>
#include <iomanip>

#ifdef ENABLE_GSL

#include <gsl/gsl_sf_dilog.h>

namespace {

double gsl_Li2(double x) {
   gsl_sf_result li2_gsl{};
   gsl_sf_dilog_e(x, &li2_gsl);
   return li2_gsl.val;
}

std::complex<double> gsl_Li2(std::complex<double> z) {
   gsl_sf_result li2_gsl_re{}, li2_gsl_im{};
   gsl_sf_complex_dilog_e(std::abs(z), std::arg(z), &li2_gsl_re, &li2_gsl_im);
   return {li2_gsl_re.val, li2_gsl_im.val};
}

} // anonymous namespace

#endif

namespace {

std::complex<long double> lt_Li2(std::complex<long double> z) {
   long double re{}, im{};
   feynhiggs_dilog(std::real(z), std::imag(z), &re, &im);
   return { re, im };
}

std::complex<double> hdecay_Li2(std::complex<double> z) {
   double li2_re{}, li2_im{};
   hdecay_dilog(std::real(z), std::imag(z), &li2_re, &li2_im);
   return { li2_re, li2_im };
}

std::complex<double> hollik_Li2(std::complex<double> z) {
   double li2_re{}, li2_im{};
   hollik_dilog(std::real(z), std::imag(z), &li2_re, &li2_im);
   return { li2_re, li2_im };
}

inline float poly_Li2(float z) {
   return li2f(z);
}

inline double poly_Li2(double z) {
   return li2(z);
}

inline long double poly_Li2(long double z) {
   return li2l(z);
}

std::complex<double> poly_Li2(std::complex<double> z) {
   double re{}, im{};
   cli2_c(std::real(z), std::imag(z), &re, &im);
   return { re, im };
}

std::complex<long double> poly_Li2(std::complex<long double> z) {
   long double re{}, im{};
   cli2l_c(std::real(z), std::imag(z), &re, &im);
   return { re, im };
}

#ifdef ENABLE_FORTRAN

double poly_Li2_fortran(double x) {
   double res{};
   li2_fortran(&x, &res);
   return res;
}

double poly_Li3_fortran(double x) {
   double res{};
   li3_fortran(&x, &res);
   return res;
}

double poly_Li4_fortran(double x) {
   double res{};
   li4_fortran(&x, &res);
   return res;
}

std::complex<double> poly_Li2_fortran(std::complex<double> z) {
   const double re = std::real(z);
   const double im = std::imag(z);
   double res_re{}, res_im{};
   cli2_fortran(&re, &im, &res_re, &res_im);
   return { res_re, res_im };
}

std::complex<double> poly_Li3_fortran(std::complex<double> z) {
   const double re = std::real(z);
   const double im = std::imag(z);
   double res_re{}, res_im{};
   cli3_fortran(&re, &im, &res_re, &res_im);
   return { res_re, res_im };
}

std::complex<double> poly_Li4_fortran(std::complex<double> z) {
   const double re = std::real(z);
   const double im = std::imag(z);
   double res_re{}, res_im{};
   cli4_fortran(&re, &im, &res_re, &res_im);
   return { res_re, res_im };
}

std::complex<double> poly_Li5_fortran(std::complex<double> z) {
   const double re = std::real(z);
   const double im = std::imag(z);
   double res_re{}, res_im{};
   cli5_fortran(&re, &im, &res_re, &res_im);
   return { res_re, res_im };
}

std::complex<double> poly_Li6_fortran(std::complex<double> z) {
   const double re = std::real(z);
   const double im = std::imag(z);
   double res_re{}, res_im{};
   cli6_fortran(&re, &im, &res_re, &res_im);
   return { res_re, res_im };
}

#endif

std::complex<double> poly_Li3(std::complex<double> z) {
   double re{}, im{};
   cli3_c(std::real(z), std::imag(z), &re, &im);
   return { re, im };
}

std::complex<long double> poly_Li3(std::complex<long double> z) {
   long double re{}, im{};
   cli3l_c(std::real(z), std::imag(z), &re, &im);
   return { re, im };
}

std::complex<double> poly_Li4(std::complex<double> z) {
   double re{}, im{};
   cli4_c(std::real(z), std::imag(z), &re, &im);
   return { re, im };
}

std::complex<long double> poly_Li4(std::complex<long double> z) {
   long double re{}, im{};
   cli4l_c(std::real(z), std::imag(z), &re, &im);
   return { re, im };
}

std::complex<double> poly_Li5(std::complex<double> z) {
   double re{}, im{};
   cli5_c(std::real(z), std::imag(z), &re, &im);
   return { re, im };
}

std::complex<long double> poly_Li5(std::complex<long double> z) {
   long double re{}, im{};
   cli5l_c(std::real(z), std::imag(z), &re, &im);
   return { re, im };
}

std::complex<double> poly_Li6(std::complex<double> z) {
   double re{}, im{};
   cli6_c(std::real(z), std::imag(z), &re, &im);
   return { re, im };
}

std::complex<long double> poly_Li6(std::complex<long double> z) {
   long double re{}, im{};
   cli6l_c(std::real(z), std::imag(z), &re, &im);
   return { re, im };
}

std::complex<double> sherpa_Li2(std::complex<double> z) {
   double re{}, im{};
   sherpa_dilog(std::real(z), std::imag(z), &re, &im);
   return { re, im };
}

std::complex<double> spheno_Li2(std::complex<double> z) {
   double re{}, im{};
   spheno_dilog(std::real(z), std::imag(z), &re, &im);
   return { re, im };
}

long double tsil_Li2(long double x)
{
   return tsil_dilog_real(x);
}

std::complex<long double> tsil_Li2(std::complex<long double> z)
{
   long double re{}, im{};
   tsil_dilog_complex(std::real(z), std::imag(z), &re, &im);
   return { re, im };
}

std::complex<long double> tsil_Li3(std::complex<long double> z)
{
   long double re{}, im{};
   tsil_trilog_complex(std::real(z), std::imag(z), &re, &im);
   return { re, im };
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
   const auto min = 0.0;
   const auto max = 0.5;

   const auto values_f  = generate_random_scalars<float>(N, min, max);
   const auto values_d  = generate_random_scalars<double>(N, min, max);
   const auto values_l  = generate_random_scalars<long double>(N, min, max);
   const auto values_cf = generate_random_complexes<float>(N, min, max);
   const auto values_cd = generate_random_complexes<double>(N, min, max);
   const auto values_cl = generate_random_complexes<long double>(N, min, max);

   print_headline("Li2 (real)");

   bench_fn([&](float x) { return polylogarithm::Li2(x); }, values_f,
            "polylogarithm C++", "float");

   bench_fn([&](float x) { return poly_Li2(x); }, values_f,
            "polylogarithm C", "float");

   bench_fn([&](double x) { return polylogarithm::Li2(x); }, values_d,
            "polylogarithm C++", "double");

   bench_fn([&](double x) { return poly_Li2(x); }, values_d,
            "polylogarithm C", "double");

#ifdef ENABLE_FORTRAN
   bench_fn([&](double x) { return poly_Li2_fortran(x); }, values_d,
            "polylogarithm Fortran", "double");
#endif

#ifdef ENABLE_GSL
   bench_fn([&](double x) { return gsl_Li2(x); }, values_d,
            "GSL", "double");
#endif

   bench_fn([&](double x) { return babar_dilog(x); }, values_d,
            "babar", "double");

   bench_fn([&](double x) { return cephes_dilog(x); }, values_d,
            "cephes", "double");

   bench_fn([&](double x) { return cephes_dilog_2(x); }, values_d,
            "cephes 2", "double");

   bench_fn([&](double x) { return algorithm_327(x); }, values_d,
            "algorithm 327", "double");

   bench_fn([&](double x) { return algorithm_490(x); }, values_d,
            "algorithm 490", "double");

   bench_fn([&](double x) { return hassani_dilog(x); }, values_d,
            "Hassani", "double");

   bench_fn([&](double x) { return koelbig_dilog(x); }, values_d,
            "Koelbig", "double");

   bench_fn([&](double x) { return morris_dilog(x); }, values_d,
            "Morris", "double");

   bench_fn([&](long double x) { return polylogarithm::Li2(x); }, values_l,
            "polylogarithm C++", "long double");

   bench_fn([&](long double x) { return poly_Li2(x); }, values_l,
            "polylogarithm C", "long double");

   bench_fn([&](long double x) { return koelbig_dilogl(x); }, values_l,
            "Koelbig", "long double");

   bench_fn([&](long double x) { return tsil_Li2(x); }, values_l,
            "TSIL", "long double");

   print_headline("Li2 (complex)");

   bench_fn([&](std::complex<float> z) { return polylogarithm::Li2(z); },
            values_cf, "polylogarithm C++", "float");

   bench_fn([&](std::complex<double> z) { return polylogarithm::Li2(z); },
            values_cd, "polylogarithm C++", "double");

   bench_fn([&](std::complex<double> z) { return poly_Li2(z); },
            values_cd, "polylogarithm C", "double");

#ifdef ENABLE_FORTRAN
   bench_fn([&](std::complex<double> z) { return poly_Li2_fortran(z); },
            values_cd, "polylogarithm Fortran", "double");
#endif

#ifdef ENABLE_GSL
   bench_fn([&](std::complex<double> z) { return gsl_Li2(z); },
            values_cd, "GSL", "double");
#endif

   bench_fn([&](std::complex<double> z) { return hdecay_Li2(z); },
            values_cd, "HDECAY", "double");

   bench_fn([&](std::complex<double> z) { return hollik_Li2(z); },
            values_cd, "Hollik", "double");

   bench_fn([&](std::complex<double> z) { return sherpa_Li2(z); },
            values_cd, "Sherpa", "double");

   bench_fn([&](std::complex<double> z) { return spheno_Li2(z); },
            values_cd, "SPheno", "double");

   bench_fn([&](std::complex<long double> z) { return polylogarithm::Li2(z); },
            values_cl, "polylogarithm C++",  "long double");

   bench_fn([&](std::complex<long double> z) { return poly_Li2(z); },
            values_cl, "polylogarithm C",  "long double");

   bench_fn([&](std::complex<long double> z) { return lt_Li2(z); },
            values_cl, "LoopTools", "long double");

   bench_fn([&](std::complex<long double> z) { return tsil_Li2(z); },
            values_cl, "TSIL", "long double");

   print_headline("Li3 (real)");

   bench_fn([&](double x) { return polylogarithm::Li3(x); }, values_d,
            "polylogarithm C++", "double");

   bench_fn([&](double x) { return li3(x); }, values_d,
            "polylogarithm C", "double");

#ifdef ENABLE_FORTRAN
   bench_fn([&](double x) { return poly_Li3_fortran(x); }, values_d,
            "polylogarithm Fortran", "double");
#endif

   print_headline("Li3 (complex)");

   bench_fn([&](std::complex<double> z) { return polylogarithm::Li3(z); },
            values_cd, "polylogarithm C++", "double");

   bench_fn([&](std::complex<double> z) { return poly_Li3(z); },
            values_cd, "polylogarithm C", "double");

#ifdef ENABLE_FORTRAN
   bench_fn([&](std::complex<double> z) { return poly_Li3_fortran(z); },
            values_cd, "polylogarithm Fortran", "double");
#endif

   bench_fn([&](std::complex<long double> z) { return polylogarithm::Li3(z); },
            values_cl, "polylogarithm C++", "long double");

   bench_fn([&](std::complex<long double> z) { return poly_Li3(z); },
            values_cl, "polylogarithm C", "long double");

   bench_fn([&](std::complex<long double> z) { return tsil_Li3(z); },
            values_cl, "TSIL", "long double");

   print_headline("Li4 (real)");

   bench_fn([&](double x) { return polylogarithm::Li4(x); }, values_d,
            "polylogarithm C++", "double");

   bench_fn([&](double x) { return li4(x); }, values_d,
            "polylogarithm C", "double");

#ifdef ENABLE_FORTRAN
   bench_fn([&](double x) { return poly_Li4_fortran(x); }, values_d,
            "polylogarithm Fortran", "double");
#endif

   print_headline("Li4 (complex)");

   bench_fn([&](std::complex<double> z) { return polylogarithm::Li4(z); },
            values_cd, "polylogarithm C++", "double");

   bench_fn([&](std::complex<double> z) { return poly_Li4(z); },
            values_cd, "polylogarithm C", "double");

#ifdef ENABLE_FORTRAN
   bench_fn([&](std::complex<double> z) { return poly_Li4_fortran(z); },
            values_cd, "polylogarithm Fortran", "double");
#endif

   bench_fn([&](std::complex<long double> z) { return polylogarithm::Li4(z); },
            values_cd, "polylogarithm C++", "long double");

   bench_fn([&](std::complex<long double> z) { return poly_Li4(z); },
            values_cd, "polylogarithm C", "long double");

   print_headline("Li5 (complex)");

   bench_fn([&](std::complex<double> z) { return polylogarithm::Li5(z); },
            values_cd, "polylogarithm C++", "double");

   bench_fn([&](std::complex<double> z) { return poly_Li5(z); },
            values_cd, "polylogarithm C", "double");

#ifdef ENABLE_FORTRAN
   bench_fn([&](std::complex<double> z) { return poly_Li5_fortran(z); },
            values_cd, "polylogarithm Fortran", "double");
#endif

   bench_fn([&](std::complex<long double> z) { return polylogarithm::Li5(z); },
            values_cd, "polylogarithm C++", "long double");

   bench_fn([&](std::complex<long double> z) { return poly_Li5(z); },
            values_cd, "polylogarithm C", "long double");

   print_headline("Li6 (complex)");

   bench_fn([&](std::complex<double> z) { return polylogarithm::Li6(z); },
            values_cd, "polylogarithm C++", "double");

   bench_fn([&](std::complex<double> z) { return poly_Li6(z); },
            values_cd, "polylogarithm C", "double");

#ifdef ENABLE_FORTRAN
   bench_fn([&](std::complex<double> z) { return poly_Li6_fortran(z); },
            values_cd, "polylogarithm Fortran", "double");
#endif

   bench_fn([&](std::complex<long double> z) { return polylogarithm::Li6(z); },
            values_cd, "polylogarithm C++", "long double");

   bench_fn([&](std::complex<long double> z) { return poly_Li6(z); },
            values_cd, "polylogarithm C", "long double");

   print_headline("Li(n,z) (complex)");

   const auto values_cd_small = generate_random_complexes<double>(N/10, min, max);

   bench_fn([&](std::complex<double> z) { return polylogarithm::Li(-1000,z); },
            values_cd_small, "n=-1000 polylogarithm", "double");

   bench_fn([&](std::complex<double> z) { return polylogarithm::Li(-100,z); },
            values_cd_small, "n=-100 polylogarithm", "double");

   bench_fn([&](std::complex<double> z) { return polylogarithm::Li(-10,z); },
            values_cd_small, "n=-10  polylogarithm", "double");

   bench_fn([&](std::complex<double> z) { return polylogarithm::Li(-6,z); },
            values_cd_small, "n=-6   polylogarithm", "double");

   bench_fn([&](std::complex<double> z) { return polylogarithm::Li(6,z); },
            values_cd_small, "n=6    polylogarithm", "double");

   bench_fn([&](std::complex<double> z) { return polylogarithm::Li(10,z); },
            values_cd_small, "n=10   polylogarithm", "double");

   bench_fn([&](std::complex<double> z) { return polylogarithm::Li(100,z); },
            values_cd_small, "n=100  polylogarithm", "double");

   bench_fn([&](std::complex<double> z) { return polylogarithm::Li(1000,z); },
            values_cd_small, "n=1000  polylogarithm", "double");

   bench_fn([&](std::complex<double> z) { return polylogarithm::Li(10000,z); },
            values_cd_small, "n=10000  polylogarithm", "double");

   return 0;
}
