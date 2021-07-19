#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "alt.h"
#include "Cl2.hpp"
#include "Li2.hpp"
#include "read_data.hpp"
#include <cmath>
#include <complex>
#include <limits>
#include <vector>

#define CHECK_CLOSE(a,b,eps) CHECK((a) == doctest::Approx(b).epsilon(eps))
#define CHECK_SMALL(a,eps) CHECK(std::abs(a) <= (eps))

#ifdef ENABLE_GSL

#include <gsl/gsl_sf_clausen.h>

double gsl_cl2(double x) {
   return gsl_sf_clausen(x);;
}

#endif

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

TEST_CASE("test_real_fixed_values")
{
   const double pi64 = 3.1415926535897932;
   const double pi128 = 3.14159265358979323846264338327950288419717;
   const auto eps64  = std::pow(10.0 , -std::numeric_limits<double>::digits10);
   const auto eps128 = std::pow(10.0L, -std::numeric_limits<long double>::digits10);
   const std::string filename(std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "Cl2.txt");
   const auto fixed_values = polylogarithm::test::read_reals_from_file<long double>(filename);

   for (auto v: fixed_values) {
      const auto x128 = v.first;
      const auto x64 = static_cast<double>(x128);
      const auto cl128_expected = v.second;
      const auto cl64_expected = static_cast<double>(cl128_expected);

      const auto cl64_koelbig = koelbig_cl2(x64);
      const auto cl64_poly    = polylogarithm::Cl2(x64);
      const auto cl128_poly   = polylogarithm::Cl2(x128);

#ifdef ENABLE_GSL
      const auto cl64_gsl     = gsl_cl2(x64);
#endif

      INFO("x(64)         = " << x64);
      INFO("Cl2(64)  real = " << cl64_expected  << " (expected)");
      INFO("Cl2(64)  real = " << cl64_poly      << " (polylogarithm C++)");
#ifdef ENABLE_GSL
      INFO("Cl2(64)  cmpl = " << cl64_gsl       << " (GSL)");
#endif
      INFO("Cl2(64)  cmpl = " << cl64_koelbig   << " (Koelbig)");
      INFO("------------------------------------------------------------");
      INFO("x(128)        = " << x128);
      INFO("Cl2(128) cmpl = " << cl128_expected << " (expected)");
      INFO("Cl2(128) cmpl = " << cl128_poly     << " (polylogarithm C++)");

      if (std::abs(x64 - 2*pi64) > 1e-2) {
         CHECK_CLOSE(cl64_poly   , cl64_expected , 2*eps64);
      } else if (std::abs(x64 - 2*pi64) > 1e-12) {
         CHECK_CLOSE(cl64_poly   , cl64_expected , 10*eps64);
      } else {
         CHECK_CLOSE(cl64_poly   , cl64_expected , 100*eps64);
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
      if (std::abs(x128 - 2*pi128) > 1e-10L) {
         CHECK_CLOSE(cl128_poly  , cl128_expected, 6*eps128);
      } else {
         CHECK_CLOSE(cl128_poly  , cl128_expected, 50*eps128);
      }
   }
}
