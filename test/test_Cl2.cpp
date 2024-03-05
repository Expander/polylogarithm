#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "alt.h"
#include "c_wrappers.h"
#include "fortran_wrappers.h"
#include "Cl2.hpp"
#include "Li2.hpp"
#include "read_data.hpp"
#include "test.hpp"
#include <cmath>
#include <complex>
#include <limits>
#include <vector>

#ifdef ENABLE_GSL
#include <gsl/gsl_sf_clausen.h>
#endif

namespace {

double poly_Cl2(double x) {
   return cl2(x);
}

long double poly_Cl2(long double x) {
   return cl2l(x);
}

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

long double Cl2_via_Li2(long double x) noexcept
{
   return std::imag(polylogarithm::Li2(std::polar(1.0L, x)));
}

} // anonymous namespace

std::vector<double> float_range(
   double start, double stop, std::size_t number_of_steps)
{
   const double step_size = (stop - start) / number_of_steps;
   std::vector<double> result(number_of_steps);

   for (std::size_t i = 0; i < number_of_steps; ++i) {
      const double point = start + i * step_size;
      result[i] = point;
   }

   return result;
}

TEST_CASE("test_special_values")
{
   using polylogarithm::Cl2;
   const double pi  = M_PI;
   const double catalan = 0.91596559417721901505460351493238411077414937428167;

   CHECK_CLOSE(Cl2(pi/2.), catalan, 1e-15);
}

TEST_CASE("test_kummer_relation")
{
   using polylogarithm::Cl2;
   using polylogarithm::Li2;
   const double pi  = M_PI;
   const double z2  = 1.644934066848226436472415166646025189218949901206798437735558229;
   const std::complex<double> i(0.,1.);

   const auto thetas = float_range(0., 2*pi, 100);

   for (const auto t: thetas) {
      const auto lhs = Li2(std::exp(i*t));
      const auto rhs = z2 - t*(2*pi - t)/4. + i*Cl2(t);

      CHECK_CLOSE(std::real(lhs), std::real(rhs), 1e-15);
      CHECK_CLOSE(std::imag(lhs), std::imag(rhs), 1e-15);
   }
}

TEST_CASE("test_roots")
{
   using polylogarithm::Cl2;
   const double pi  = M_PI;

   for (int k = -10; k < 10; k++) {
      CHECK_SMALL(Cl2(k*pi), 1e-10);
   }
}

// tests signbit for 0.0 and -0.0 arguments
TEST_CASE("test_signed_zero")
{
   // skip test if platform does not supprt signed zero
   if (!has_signed_zero()) {
      return;
   }

   using polylogarithm::Cl2;

   const float  pz32 = 0.0f, nz32 = -0.0f;
   const double pz64 = 0.0, nz64 = -0.0;
   const long double pz128 = 0.0L, nz128 = -0.0L;

   CHECK( std::signbit(Cl2(nz32)));
   CHECK(!std::signbit(Cl2(pz32)));
   CHECK( std::signbit(poly_Cl2(nz32)));
   CHECK(!std::signbit(poly_Cl2(pz32)));

   CHECK( std::signbit(Cl2(nz64)));
   CHECK(!std::signbit(Cl2(pz64)));
   CHECK( std::signbit(poly_Cl2(nz64)));
   CHECK(!std::signbit(poly_Cl2(pz64)));
#ifdef ENABLE_FORTRAN
   CHECK( std::signbit(poly_Cl2_fortran(nz64)));
   CHECK(!std::signbit(poly_Cl2_fortran(pz64)));
#endif

   CHECK( std::signbit(Cl2(nz128)));
   CHECK(!std::signbit(Cl2(pz128)));
   CHECK( std::signbit(poly_Cl2(nz128)));
   CHECK(!std::signbit(poly_Cl2(pz128)));
}

TEST_CASE("test_real_fixed_values")
{
   const double pi64 = 3.1415926535897932;
   const long double pi128 = 3.14159265358979323846264338327950288419717L;
   const auto eps64  = std::pow(10.0 , -std::numeric_limits<double>::digits10);
   const auto eps128 = std::pow(10.0L, -std::numeric_limits<long double>::digits10);
   const std::string filename(std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "Cl2.txt");
   const auto fixed_values = polylogarithm::test::read_reals_from_file<long double>(filename);

   for (auto v: fixed_values) {
      const auto x128 = v.first;
      const auto x64 = static_cast<double>(x128);
      const auto cl128_expected = v.second;
      const auto cl64_expected = static_cast<double>(cl128_expected);

      const auto cl64_bernoulli = clausen_2_bernoulli(x64);
      const auto cl64_koelbig   = clausen_2_koelbig(x64);
      const auto cl64_pade      = clausen_2_pade(x64);
      const auto cl64_poly      = polylogarithm::Cl2(x64);
      const auto cl64_poly_c    = cl2(x64);
#ifdef ENABLE_FORTRAN
      const auto cl64_poly_f    = poly_Cl2_fortran(x64);
#endif
      const auto cl64_li2       = Cl2_via_Li2(x64);
      const auto cl64_wu        = clausen_2_wu(x64);
      const auto cl128_poly     = polylogarithm::Cl2(x128);
      const auto cl128_poly_c   = cl2l(x128);
      const auto cl128_koelbig  = clausen_2l_koelbig(x128);
      const auto cl128_li2      = Cl2_via_Li2(x128);

#ifdef ENABLE_GSL
      const auto cl64_gsl     = gsl_sf_clausen(x64);
#endif

      INFO("x(64)         = " << x64);
      INFO("Cl2(64)  real = " << cl64_expected  << " (expected)");
      INFO("Cl2(64)  real = " << cl64_poly      << " (polylogarithm C++)");
      INFO("Cl2(64)  real = " << cl64_poly_c    << " (polylogarithm C)");
#ifdef ENABLE_FORTRAN
      INFO("Cl2(64)  real = " << cl64_poly_f    << " (polylogarithm Fortran)");
#endif
      INFO("Cl2(64)  real = " << cl64_li2       << " (via Li2 C++)");
#ifdef ENABLE_GSL
      INFO("Cl2(64)  real = " << cl64_gsl       << " (GSL)");
#endif
      INFO("Cl2(64)  real = " << cl64_bernoulli << " (Bernoulli series)");
      INFO("Cl2(64)  real = " << cl64_koelbig   << " (Koelbig)");
      INFO("Cl2(64)  real = " << cl64_pade      << " (Pade)");
      INFO("Cl2(64)  real = " << cl64_wu        << " (Wu et.al.)");
      INFO("------------------------------------------------------------");
      INFO("x(128)        = " << x128);
      INFO("Cl2(128) real = " << cl128_expected << " (expected)");
      INFO("Cl2(128) real = " << cl128_poly     << " (polylogarithm C++)");
      INFO("Cl2(128) real = " << cl128_poly_c   << " (polylogarithm C)");
      INFO("Cl2(128) real = " << cl128_koelbig  << " (Koelbig)");
      INFO("Cl2(128) real = " << cl128_li2      << " (via Li2 C++)");

      if (std::abs(x64 - 2*pi64) > 1e-2) {
         CHECK_CLOSE(cl64_poly   , cl64_expected , 2*eps64);
      } else if (std::abs(x64 - 2*pi64) > 1e-12) {
         CHECK_CLOSE(cl64_poly   , cl64_expected , 10*eps64);
      } else {
         CHECK_CLOSE(cl64_poly   , cl64_expected , 100*eps64);
      }
      if (std::abs(x64 - 2*pi64) > 1e-2) {
         CHECK_CLOSE(cl64_poly_c , cl64_expected , 2*eps64);
      } else if (std::abs(x64 - 2*pi64) > 1e-12) {
         CHECK_CLOSE(cl64_poly_c , cl64_expected , 10*eps64);
      } else {
         CHECK_CLOSE(cl64_poly_c , cl64_expected , 100*eps64);
      }
#ifdef ENABLE_FORTRAN
      if (std::abs(x64 - 2*pi64) > 1e-2) {
         CHECK_CLOSE(cl64_poly_f , cl64_expected , 2*eps64);
      } else if (std::abs(x64 - 2*pi64) > 1e-12) {
         CHECK_CLOSE(cl64_poly_f , cl64_expected , 10*eps64);
      } else {
         CHECK_CLOSE(cl64_poly_f , cl64_expected , 100*eps64);
      }
#endif
      CHECK_CLOSE(cl64_li2       , cl64_expected , 10*eps64);
      if (std::abs(x64 - 2*pi64) > 1e-3) {
         CHECK_CLOSE(cl64_bernoulli, cl64_expected , 2*eps64);
      } else {
         CHECK_CLOSE(cl64_bernoulli, cl64_expected , 10*eps64);
      }
#ifdef ENABLE_GSL
      if (std::abs(x64 - 2*pi64) > 1e-3) {
         CHECK_CLOSE(cl64_gsl    , cl64_expected , 2*eps64);
      } else {
         CHECK_CLOSE(cl64_gsl    , cl64_expected , 10*eps64);
      }
#endif
      if (std::abs(x64 - 2*pi64) > 1e-2) {
         CHECK_CLOSE(cl64_koelbig, cl64_expected , 2*eps64);
      } else if (std::abs(x64 - 2*pi64) > 1e-12) {
         CHECK_CLOSE(cl64_koelbig, cl64_expected , 10*eps64);
      } else {
         CHECK_CLOSE(cl64_koelbig, cl64_expected , 100*eps64);
      }
      if (std::abs(x64 - 2*pi64) > 1e-2) {
         CHECK_CLOSE(cl64_pade, cl64_expected , 2*eps64);
      } else if (std::abs(x64 - 2*pi64) > 1e-12) {
         CHECK_CLOSE(cl64_pade, cl64_expected , 10*eps64);
      } else {
         CHECK_CLOSE(cl64_pade, cl64_expected , 100*eps64);
      }
      if (std::abs(x64 - 2*pi64) > 1e-3) {
         CHECK_CLOSE(cl64_wu  , cl64_expected , 2*eps64);
      } else {
         CHECK_CLOSE(cl64_wu  , cl64_expected , 10*eps64);
      }

      if (std::abs(x128 - 2*pi128) > 1e-10L) {
         CHECK_CLOSE(cl128_poly  , cl128_expected, 6*eps128);
      } else {
         CHECK_CLOSE(cl128_poly  , cl128_expected, 50*eps128);
      }
      if (std::abs(x128 - 2*pi128) > 1e-10L) {
         CHECK_CLOSE(cl128_poly_c, cl128_expected, 6*eps128);
      } else {
         CHECK_CLOSE(cl128_poly_c, cl128_expected, 50*eps128);
      }
      if (std::abs(x128 - 2*pi128) > 1e-10L) {
         CHECK_CLOSE(cl128_koelbig, cl128_expected, 6*eps128);
      } else {
         CHECK_CLOSE(cl128_koelbig, cl128_expected, 50*eps128);
      }
      CHECK_CLOSE(cl128_li2      , cl128_expected, 10*eps128);
   }
}
