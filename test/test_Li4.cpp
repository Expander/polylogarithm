#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "alt.h"
#include "c_wrappers.h"
#include "fortran_wrappers.h"
#include "Li4.hpp"
#include "read_data.hpp"
#include "test.hpp"
#include <cmath>
#include <limits>
#include <utility>

#ifndef M_PI
#define M_PI 3.1415926535897932
#endif

#define CHECK_CLOSE(a,b,eps) CHECK((a) == doctest::Approx(b).epsilon(eps))
#define CHECK_CLOSE_COMPLEX(a,b,eps) do {                               \
      CHECK_CLOSE(std::real(a), std::real(b), (eps));                   \
      CHECK_CLOSE(std::imag(a), std::imag(b), (eps));                   \
   } while (0)
#define CHECK_SMALL(a,eps) CHECK(std::abs(a) < (eps))

template <class T> T sqr(T x) { return x*x; }
template <class T> T pow4(T x) { return x*x*x*x; }

double poly_Li4(double z) {
   return li4(z);
}

std::complex<double> poly_Li4(std::complex<double> z) {
   double re{}, im{};
   cli4_c(std::real(z), std::imag(z), &re, &im);
   return { re, im };
}

#ifdef ENABLE_FORTRAN

double poly_Li4_fortran(double x) {
   double res{};
   li4_fortran(&x, &res);
   return res;
}

std::complex<double> poly_Li4_fortran(std::complex<double> z) {
   const double re = std::real(z);
   const double im = std::imag(z);
   double res_re{}, res_im{};
   cli4_fortran(&re, &im, &res_re, &res_im);
   return { res_re, res_im };
}

#endif

std::complex<long double> poly_Li4(std::complex<long double> z) {
   long double re{}, im{};
   cli4l_c(std::real(z), std::imag(z), &re, &im);
   return { re, im };
}

TEST_CASE("test_special_values")
{
   using polylogarithm::Li4;

   const double eps = std::pow(10.0, -std::numeric_limits<double>::digits10);
   const double pi = M_PI;
   const double zeta4 = 1.082323233711138;
   const std::complex<double> zero(0.0, 0.0);
   const std::complex<double> one(1.0, 0.0);
   const std::complex<double> mone(-1.0, 0.0);
   const std::complex<double> half(0.5, 0.0);

   CHECK_CLOSE_COMPLEX(Li4(zero), 0, eps);
   CHECK_CLOSE_COMPLEX(Li4(one), zeta4, eps);
   CHECK_CLOSE_COMPLEX(Li4(mone), -7.*pow4(pi)/720, eps);
   CHECK_CLOSE_COMPLEX(Li4(half), 0.5174790616738994, eps);
}

// tests signbit for 0.0 and -0.0 arguments
TEST_CASE("test_signed_zero")
{
   // skip test if platform does not supprt signed zero
   if (!has_signed_zero()) {
      return;
   }

   using polylogarithm::Li4;

   const float  pz32 = 0.0f, nz32 = -0.0f;
   const double pz64 = 0.0, nz64 = -0.0;
   const long double pz128 = 0.0L, nz128 = -0.0L;

   // real Li4
   CHECK( std::signbit(Li4(nz64)));
   CHECK(!std::signbit(Li4(pz64)));
   CHECK( std::signbit(poly_Li4(nz64)));
   CHECK(!std::signbit(poly_Li4(pz64)));
#ifdef ENABLE_FORTRAN
   CHECK( std::signbit(poly_Li4_fortran(nz64)));
   CHECK(!std::signbit(poly_Li4_fortran(pz64)));
#endif

   // complex Li4
   CHECK( std::signbit(std::real(Li4(std::complex<double>(nz64, nz64)))));
   CHECK( std::signbit(std::imag(Li4(std::complex<double>(nz64, nz64)))));
   CHECK(!std::signbit(std::real(Li4(std::complex<double>(pz64, nz64)))));
   CHECK( std::signbit(std::imag(Li4(std::complex<double>(pz64, nz64)))));
   CHECK( std::signbit(std::real(Li4(std::complex<double>(nz64, pz64)))));
   CHECK(!std::signbit(std::imag(Li4(std::complex<double>(nz64, pz64)))));
   CHECK(!std::signbit(std::real(Li4(std::complex<double>(pz64, pz64)))));
   CHECK(!std::signbit(std::imag(Li4(std::complex<double>(pz64, pz64)))));

   CHECK( std::signbit(std::real(poly_Li4(std::complex<double>(nz64, nz64)))));
   CHECK( std::signbit(std::imag(poly_Li4(std::complex<double>(nz64, nz64)))));
   CHECK(!std::signbit(std::real(poly_Li4(std::complex<double>(pz64, nz64)))));
   CHECK( std::signbit(std::imag(poly_Li4(std::complex<double>(pz64, nz64)))));
   CHECK( std::signbit(std::real(poly_Li4(std::complex<double>(nz64, pz64)))));
   CHECK(!std::signbit(std::imag(poly_Li4(std::complex<double>(nz64, pz64)))));
   CHECK(!std::signbit(std::real(poly_Li4(std::complex<double>(pz64, pz64)))));
   CHECK(!std::signbit(std::imag(poly_Li4(std::complex<double>(pz64, pz64)))));

#ifdef ENABLE_FORTRAN
   CHECK( std::signbit(std::real(poly_Li4_fortran(std::complex<double>(nz64, nz64)))));
   CHECK( std::signbit(std::imag(poly_Li4_fortran(std::complex<double>(nz64, nz64)))));
   CHECK(!std::signbit(std::real(poly_Li4_fortran(std::complex<double>(pz64, nz64)))));
   CHECK( std::signbit(std::imag(poly_Li4_fortran(std::complex<double>(pz64, nz64)))));
   CHECK( std::signbit(std::real(poly_Li4_fortran(std::complex<double>(nz64, pz64)))));
   CHECK(!std::signbit(std::imag(poly_Li4_fortran(std::complex<double>(nz64, pz64)))));
   CHECK(!std::signbit(std::real(poly_Li4_fortran(std::complex<double>(pz64, pz64)))));
   CHECK(!std::signbit(std::imag(poly_Li4_fortran(std::complex<double>(pz64, pz64)))));
#endif

   CHECK( std::signbit(std::real(Li4(std::complex<long double>(nz128, nz128)))));
   CHECK( std::signbit(std::imag(Li4(std::complex<long double>(nz128, nz128)))));
   CHECK(!std::signbit(std::real(Li4(std::complex<long double>(pz128, nz128)))));
   CHECK( std::signbit(std::imag(Li4(std::complex<long double>(pz128, nz128)))));
   CHECK( std::signbit(std::real(Li4(std::complex<long double>(nz128, pz128)))));
   CHECK(!std::signbit(std::imag(Li4(std::complex<long double>(nz128, pz128)))));
   CHECK(!std::signbit(std::real(Li4(std::complex<long double>(pz128, pz128)))));
   CHECK(!std::signbit(std::imag(Li4(std::complex<long double>(pz128, pz128)))));

   CHECK( std::signbit(std::real(poly_Li4(std::complex<long double>(nz128, nz128)))));
   CHECK( std::signbit(std::imag(poly_Li4(std::complex<long double>(nz128, nz128)))));
   CHECK(!std::signbit(std::real(poly_Li4(std::complex<long double>(pz128, nz128)))));
   CHECK( std::signbit(std::imag(poly_Li4(std::complex<long double>(pz128, nz128)))));
   CHECK( std::signbit(std::real(poly_Li4(std::complex<long double>(nz128, pz128)))));
   CHECK(!std::signbit(std::imag(poly_Li4(std::complex<long double>(nz128, pz128)))));
   CHECK(!std::signbit(std::real(poly_Li4(std::complex<long double>(pz128, pz128)))));
   CHECK(!std::signbit(std::imag(poly_Li4(std::complex<long double>(pz128, pz128)))));
}

template<typename T>
struct Data {
   std::complex<T> z;
   std::complex<T> li_expected;
   T eps;
};

TEST_CASE("test_overflow")
{
   using polylogarithm::Li4;

   const double eps64 = std::pow(10.0, -std::numeric_limits<double>::digits10);
   const long double eps128 = std::pow(10.0L, -std::numeric_limits<long double>::digits10);

   const Data<double> data64[] = {
      {{1e300, 1.0}, {-9.4863817894708364e9, 1.725875455850714e8}, eps64},
      {{1.0, 1e300}, {-9.4872648206269765e9, 8.62951114411071e7}, eps64}
   };

   for (const auto& d : data64) {
      CHECK_CLOSE_COMPLEX(Li4(d.z), d.li_expected, d.eps);
      CHECK_CLOSE_COMPLEX(poly_Li4(d.z), d.li_expected, d.eps);
#ifdef ENABLE_FORTRAN
      CHECK_CLOSE_COMPLEX(poly_Li4_fortran(d.z), d.li_expected, d.eps);
#endif
   }

   const Data<long double> data128[] = {
      {{1e4000L, 1.0L}, {-2.998411785811730437806291599010045237778e14L, 4.090964043497988426223363565328465514e11L}, eps128}
   };

   for (const auto& d : data128) {
      CHECK_CLOSE_COMPLEX(Li4(d.z), d.li_expected, d.eps);
      CHECK_CLOSE_COMPLEX(poly_Li4(d.z), d.li_expected, d.eps);
   }
}

TEST_CASE("test_real_fixed_values")
{
   const auto eps64 = std::pow(10.0 , -std::numeric_limits<double>::digits10);

   const std::string filename(std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "Li4.txt");
   const auto fixed_values = polylogarithm::test::read_from_file<long double>(filename);

   for (auto v: fixed_values) {
      const auto z128 = v.first;
      const auto z64 = to<double>(z128);
      const auto x64 = std::real(z64);
      const auto li128_expected = std::real(v.second);
      const auto li64_expected = static_cast<double>(li128_expected);

      if (std::imag(z128) == 0.0L) {
         const auto li64_poly   = polylogarithm::Li4(x64);
         const auto li64_poly_c = li4(x64);
#ifdef ENABLE_FORTRAN
         const auto li64_poly_f = poly_Li4_fortran(x64);
#endif

         INFO("x(64)         = " << x64);
         INFO("Li4(64)  real = " << li64_expected  << " (expected)");
         INFO("Li4(64)  real = " << li64_poly      << " (polylogarithm C++)");
         INFO("Li4(64)  real = " << li64_poly_c    << " (polylogarithm C)");
#ifdef ENABLE_FORTRAN
         INFO("Li4(64)  real = " << li64_poly_f    << " (polylogarithm Fortran)");
#endif

         CHECK_CLOSE(li64_poly  , li64_expected, 5*eps64);
         CHECK_CLOSE(li64_poly_c, li64_expected, 5*eps64);
#ifdef ENABLE_FORTRAN
         CHECK_CLOSE(li64_poly_f, li64_expected, 5*eps64);
#endif
      }
   }
}

TEST_CASE("test_complex_fixed_values")
{
   const auto eps64  = std::pow(10.0 , -std::numeric_limits<double>::digits10);
   const auto eps128 = std::pow(10.0L, -std::numeric_limits<long double>::digits10);

   const std::string filename(std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "Li4.txt");
   const auto fixed_values = polylogarithm::test::read_from_file<long double>(filename);

   for (auto v: fixed_values) {
      const auto z128 = v.first;
      const auto z64 = to<double>(z128);
      const auto li128_expected = v.second;
      const auto li64_expected = to<double>(li128_expected);

      const auto li64_cmpl    = polylogarithm::Li4(z64);
      const auto li128_cmpl   = polylogarithm::Li4(z128);
      const auto li64_cmpl_c  = poly_Li4(z64);
#ifdef ENABLE_FORTRAN
      const auto li64_cmpl_f  = poly_Li4_fortran(z64);
#endif
      const auto li128_cmpl_c = poly_Li4(z128);

      INFO("z(128)        = " << z128);
      INFO("Li4(64)  cmpl = " << li64_expected  << " (expected)");
      INFO("Li4(64)  cmpl = " << li64_cmpl      << " (polylogarithm C++)");
      INFO("Li4(64)  cmpl = " << li64_cmpl_c    << " (polylogarithm C)");
#ifdef ENABLE_FORTRAN
      INFO("Li4(64)  cmpl = " << li64_cmpl_f    << " (polylogarithm Fortran)");
#endif
      INFO("Li4(128) cmpl = " << li128_expected << " (expected)");
      INFO("Li4(128) cmpl = " << li128_cmpl     << " (polylogarithm C++)");
      INFO("Li4(128) cmpl = " << li128_cmpl_c   << " (polylogarithm C)");

      CHECK_CLOSE_COMPLEX(li64_cmpl   , li64_expected , 5*eps64);
      CHECK_CLOSE_COMPLEX(li64_cmpl_c , li64_expected , 5*eps64);
#ifdef ENABLE_FORTRAN
      CHECK_CLOSE_COMPLEX(li64_cmpl_f , li64_expected , 5*eps64);
#endif
      CHECK_CLOSE_COMPLEX(li128_cmpl  , li128_expected, 2*eps128);
      CHECK_CLOSE_COMPLEX(li128_cmpl_c, li128_expected, 2*eps128);
   }
}
