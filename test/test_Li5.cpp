#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "alt.h"
#include "c_wrappers.h"
#include "fortran_wrappers.h"
#include "Li5.hpp"
#include "read_data.hpp"
#include "test.hpp"
#include <cmath>
#include <limits>
#include <utility>

std::complex<double> poly_Li5(std::complex<double> z) {
   double re{}, im{};
   cli5_c(std::real(z), std::imag(z), &re, &im);
   return { re, im };
}

#ifdef ENABLE_FORTRAN

std::complex<double> poly_Li5_fortran(std::complex<double> z) {
   const double re = std::real(z);
   const double im = std::imag(z);
   double res_re{}, res_im{};
   cli5_fortran(&re, &im, &res_re, &res_im);
   return { res_re, res_im };
}

#endif

std::complex<long double> poly_Li5(std::complex<long double> z) {
   long double re{}, im{};
   cli5l_c(std::real(z), std::imag(z), &re, &im);
   return { re, im };
}

TEST_CASE("test_special_values")
{
   using polylogarithm::Li5;

   const double eps = std::pow(10.0, -std::numeric_limits<double>::digits10);
   const double zeta5 = 1.0369277551433699;
   const std::complex<double> zero(0.0, 0.0);
   const std::complex<double> one(1.0, 0.0);
   const std::complex<double> mone(-1.0, 0.0);
   const std::complex<double> half(0.5, 0.0);

   CHECK_CLOSE_COMPLEX(Li5(zero), 0, eps);
   CHECK_CLOSE_COMPLEX(Li5(one), zeta5, eps);
   CHECK_CLOSE_COMPLEX(Li5(mone), -15.*zeta5/16., eps);
   CHECK_CLOSE_COMPLEX(Li5(half), 0.5084005792422687, eps);
}

// tests signbit for 0.0 and -0.0 arguments
TEST_CASE("test_signed_zero")
{
   // skip test if platform does not supprt signed zero
   if (!has_signed_zero()) {
      return;
   }

   using polylogarithm::Li5;

   const double pz64 = 0.0, nz64 = -0.0;
   const long double pz128 = 0.0L, nz128 = -0.0L;

   // complex Li5
   CHECK( std::signbit(std::real(Li5(std::complex<double>(nz64, nz64)))));
   CHECK( std::signbit(std::imag(Li5(std::complex<double>(nz64, nz64)))));
   CHECK(!std::signbit(std::real(Li5(std::complex<double>(pz64, nz64)))));
   CHECK( std::signbit(std::imag(Li5(std::complex<double>(pz64, nz64)))));
   CHECK( std::signbit(std::real(Li5(std::complex<double>(nz64, pz64)))));
   CHECK(!std::signbit(std::imag(Li5(std::complex<double>(nz64, pz64)))));
   CHECK(!std::signbit(std::real(Li5(std::complex<double>(pz64, pz64)))));
   CHECK(!std::signbit(std::imag(Li5(std::complex<double>(pz64, pz64)))));

   CHECK( std::signbit(std::real(poly_Li5(std::complex<double>(nz64, nz64)))));
   CHECK( std::signbit(std::imag(poly_Li5(std::complex<double>(nz64, nz64)))));
   CHECK(!std::signbit(std::real(poly_Li5(std::complex<double>(pz64, nz64)))));
   CHECK( std::signbit(std::imag(poly_Li5(std::complex<double>(pz64, nz64)))));
   CHECK( std::signbit(std::real(poly_Li5(std::complex<double>(nz64, pz64)))));
   CHECK(!std::signbit(std::imag(poly_Li5(std::complex<double>(nz64, pz64)))));
   CHECK(!std::signbit(std::real(poly_Li5(std::complex<double>(pz64, pz64)))));
   CHECK(!std::signbit(std::imag(poly_Li5(std::complex<double>(pz64, pz64)))));

#ifdef ENABLE_FORTRAN
   CHECK( std::signbit(std::real(poly_Li5_fortran(std::complex<double>(nz64, nz64)))));
   CHECK( std::signbit(std::imag(poly_Li5_fortran(std::complex<double>(nz64, nz64)))));
   CHECK(!std::signbit(std::real(poly_Li5_fortran(std::complex<double>(pz64, nz64)))));
   CHECK( std::signbit(std::imag(poly_Li5_fortran(std::complex<double>(pz64, nz64)))));
   CHECK( std::signbit(std::real(poly_Li5_fortran(std::complex<double>(nz64, pz64)))));
   CHECK(!std::signbit(std::imag(poly_Li5_fortran(std::complex<double>(nz64, pz64)))));
   CHECK(!std::signbit(std::real(poly_Li5_fortran(std::complex<double>(pz64, pz64)))));
   CHECK(!std::signbit(std::imag(poly_Li5_fortran(std::complex<double>(pz64, pz64)))));
#endif

   CHECK( std::signbit(std::real(Li5(std::complex<long double>(nz128, nz128)))));
   CHECK( std::signbit(std::imag(Li5(std::complex<long double>(nz128, nz128)))));
   CHECK(!std::signbit(std::real(Li5(std::complex<long double>(pz128, nz128)))));
   CHECK( std::signbit(std::imag(Li5(std::complex<long double>(pz128, nz128)))));
   CHECK( std::signbit(std::real(Li5(std::complex<long double>(nz128, pz128)))));
   CHECK(!std::signbit(std::imag(Li5(std::complex<long double>(nz128, pz128)))));
   CHECK(!std::signbit(std::real(Li5(std::complex<long double>(pz128, pz128)))));
   CHECK(!std::signbit(std::imag(Li5(std::complex<long double>(pz128, pz128)))));

   CHECK( std::signbit(std::real(poly_Li5(std::complex<long double>(nz128, nz128)))));
   CHECK( std::signbit(std::imag(poly_Li5(std::complex<long double>(nz128, nz128)))));
   CHECK(!std::signbit(std::real(poly_Li5(std::complex<long double>(pz128, nz128)))));
   CHECK( std::signbit(std::imag(poly_Li5(std::complex<long double>(pz128, nz128)))));
   CHECK( std::signbit(std::real(poly_Li5(std::complex<long double>(nz128, pz128)))));
   CHECK(!std::signbit(std::imag(poly_Li5(std::complex<long double>(nz128, pz128)))));
   CHECK(!std::signbit(std::real(poly_Li5(std::complex<long double>(pz128, pz128)))));
   CHECK(!std::signbit(std::imag(poly_Li5(std::complex<long double>(pz128, pz128)))));
}

template<typename T>
struct Data {
   std::complex<T> z;
   std::complex<T> li_expected;
   T eps;
};

TEST_CASE("test_overflow")
{
   using polylogarithm::Li5;

   const double eps64 = std::pow(10.0, -std::numeric_limits<double>::digits10);
   const long double eps128 = std::pow(10.0L, -std::numeric_limits<long double>::digits10);

   const Data<double> data64[] = {
      {{1e300, 1.0}, {-1.3105197831948743e12, 2.980481322754618e10}, eps64},
      {{1.0, 1e300}, {-1.31072310968392418e12, 1.490286896860219e10}, eps64}
   };

   for (const auto& d : data64) {
      CHECK_CLOSE_COMPLEX(Li5(d.z), d.li_expected, d.eps);
      CHECK_CLOSE_COMPLEX(poly_Li5(d.z), d.li_expected, d.eps);
#ifdef ENABLE_FORTRAN
      CHECK_CLOSE_COMPLEX(poly_Li5_fortran(d.z), d.li_expected, d.eps);
#endif
   }

   const Data<long double> data128[] = {
      {{1e4000L, 1.0L}, {-5.523276910915025884712078063195501148879e17L, 9.419792822533112833122424859594537172e14L}, eps128}
   };

   for (const auto& d : data128) {
      CHECK_CLOSE_COMPLEX(Li5(d.z), d.li_expected, d.eps);
      CHECK_CLOSE_COMPLEX(poly_Li5(d.z), d.li_expected, d.eps);
   }
}

TEST_CASE("test_fixed_values")
{
   const auto eps64  = std::pow(10.0 , -std::numeric_limits<double>::digits10);
   const auto eps128 = std::pow(10.0L, -std::numeric_limits<long double>::digits10);

   const std::string filename(std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "Li5.txt");
   const auto fixed_values = polylogarithm::test::read_from_file<long double>(filename);

   for (auto v: fixed_values) {
      const auto z128 = v.first;
      const auto z64 = to<double>(z128);
      const auto li128_expected = v.second;
      const auto li64_expected = to<double>(li128_expected);

      const auto li64_cmpl  = polylogarithm::Li5(z64);
      const auto li128_cmpl = polylogarithm::Li5(z128);
      const auto li64_cmpl_c  = poly_Li5(z64);
#ifdef ENABLE_FORTRAN
      const auto li64_cmpl_f  = poly_Li5_fortran(z64);
#endif
      const auto li128_cmpl_c = poly_Li5(z128);

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

      CHECK_CLOSE_COMPLEX(li64_cmpl   , li64_expected , 2*eps64);
      CHECK_CLOSE_COMPLEX(li64_cmpl_c , li64_expected , 2*eps64);
#ifdef ENABLE_FORTRAN
      CHECK_CLOSE_COMPLEX(li64_cmpl_f , li64_expected , 2*eps64);
#endif
      CHECK_CLOSE_COMPLEX(li128_cmpl  , li128_expected, 2*eps128);
      CHECK_CLOSE_COMPLEX(li128_cmpl_c, li128_expected, 2*eps128);
   }
}
