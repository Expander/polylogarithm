#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "c_wrappers.h"
#include "Cl5.hpp"
#include "Li5.hpp"
#include "read_data.hpp"
#include <cmath>
#include <complex>
#include <vector>

#ifndef M_PI
#define M_PI 3.1415926535897932
#endif

#define CHECK_CLOSE(a,b,eps) CHECK((a) == doctest::Approx(b).epsilon(eps))
#define CHECK_SMALL(a,eps) CHECK(std::abs(a) <= (eps))

double Cl5_via_Li5(double x) noexcept
{
   return std::real(polylogarithm::Li5(std::polar(1.0, x)));
}

long double Cl5_via_Li5(long double x) noexcept
{
   return std::real(polylogarithm::Li5(std::polar(1.0L, x)));
}

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
   using polylogarithm::Cl5;
   const double pi = M_PI;
   const double z5 = 1.036927755143370;

   // Cl_5(Pi/2) = -((15 Zeta[5])/512)
   CHECK_CLOSE(Cl5(pi/2.), -15*z5/512., 1e-15);
   CHECK_CLOSE(Cl5(0.), z5, 1e-15);
   CHECK_CLOSE(Cl5(2*pi), z5, 1e-15);
}

TEST_CASE("test_duplication_formula")
{
   using polylogarithm::Cl5;
   const double pi  = M_PI;

   const auto thetas = float_range(0., pi, 100);

   for (const auto t: thetas) {
      const auto rel = Cl5(2*t) - (16*Cl5(t) + 16*Cl5(pi - t));
      CHECK_SMALL(rel, 1e-8);
   }
}

TEST_CASE("test_real_fixed_values")
{
   const auto eps64  = std::pow(10.0 , -std::numeric_limits<double>::digits10);
   const auto eps128 = std::pow(10.0L, -std::numeric_limits<long double>::digits10);
   const std::string filename(std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "Cl5.txt");
   const auto fixed_values = polylogarithm::test::read_reals_from_file<long double>(filename);

   for (auto v: fixed_values) {
      const auto x128 = v.first;
      const auto x64 = static_cast<double>(x128);
      const auto cl128_expected = v.second;
      const auto cl64_expected = static_cast<double>(cl128_expected);

      const auto cl64_li      = Cl5_via_Li5(x64);
      const auto cl64_poly    = polylogarithm::Cl5(x64);
      const auto cl64_poly_c  = cl5(x64);
      const auto cl128_li     = Cl5_via_Li5(x128);
      const auto cl128_poly   = polylogarithm::Cl5(x128);
      const auto cl128_poly_c = cl5l(x128);

      INFO("x(64)         = " << x64);
      INFO("Cl5(64)  real = " << cl64_expected  << " (expected)");
      INFO("Cl5(64)  real = " << cl64_poly      << " (polylogarithm C++)");
      INFO("Cl5(64)  real = " << cl64_poly_c    << " (polylogarithm C)");
      INFO("Cl5(64)  real = " << cl64_li        << " (via Li5 C++)");
      INFO("------------------------------------------------------------");
      INFO("x(128)        = " << x128);
      INFO("Cl5(128) real = " << cl128_expected << " (expected)");
      INFO("Cl5(128) real = " << cl128_poly     << " (polylogarithm C++)");
      INFO("Cl5(128) real = " << cl128_poly_c   << " (polylogarithm C)");
      INFO("Cl5(128) real = " << cl128_li       << " (via Li5 C++)");

      CHECK_CLOSE(cl64_li     , cl64_expected , 2*eps64 );
      CHECK_CLOSE(cl64_poly   , cl64_expected , 2*eps64 );
      CHECK_CLOSE(cl64_poly_c , cl64_expected , 2*eps64 );
      CHECK_CLOSE(cl128_li    , cl128_expected, 4*eps128);
      CHECK_CLOSE(cl128_poly  , cl128_expected, 7*eps128);
      CHECK_CLOSE(cl128_poly_c, cl128_expected, 7*eps128);
   }
}
