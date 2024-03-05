#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "c_wrappers.h"
#include "Cl6.hpp"
#include "fortran_wrappers.h"
#include "Li6.hpp"
#include "read_data.hpp"
#include "test.hpp"
#include <cmath>
#include <complex>
#include <vector>

#ifndef M_PI
#define M_PI 3.1415926535897932
#endif

#define CHECK_CLOSE(a,b,eps) CHECK((a) == doctest::Approx(b).epsilon(eps))
#define CHECK_SMALL(a,eps) CHECK(std::abs(a) <= (eps))

double poly_Cl6(double x) {
   return cl6(x);
}

long double poly_Cl6(long double x) {
   return cl6l(x);
}

#ifdef ENABLE_FORTRAN

double poly_Cl6_fortran(double x) {
   double res{};
   cl6_fortran(&x, &res);
   return res;
}

#endif

double Cl6_via_Li6(double x) noexcept
{
   return std::imag(polylogarithm::Li6(std::polar(1.0, x)));
}

long double Cl6_via_Li6(long double x) noexcept
{
   return std::imag(polylogarithm::Li6(std::polar(1.0L, x)));
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
   using polylogarithm::Cl6;
   const double pi  = M_PI;

   // Cl_6(Pi/2) = DirichletBeta(6)
   CHECK_CLOSE(Cl6(pi/2.), 0.998685222218438135441600787860206549678364546126514411, 1e-15);
}

TEST_CASE("test_duplication_formula")
{
   using polylogarithm::Cl6;
   const double pi  = M_PI;

   const auto thetas = float_range(0., pi, 100);

   for (const auto t: thetas) {
      const auto rel = Cl6(2*t) - (32*Cl6(t) - 32*Cl6(pi - t));
      CHECK_SMALL(rel, 1e-7);
   }
}

TEST_CASE("test_roots")
{
   using polylogarithm::Cl6;
   const double pi  = M_PI;

   for (int k = -10; k < 10; k++) {
      CHECK_SMALL(Cl6(k*pi), 1e-10);
   }
}

// tests signbit for 0.0 and -0.0 arguments
TEST_CASE("test_signed_zero")
{
   // skip test if platform does not supprt signed zero
   if (!has_signed_zero()) {
      return;
   }

   using polylogarithm::Cl6;

   const float  pz32 = 0.0f, nz32 = -0.0f;
   const double pz64 = 0.0, nz64 = -0.0;
   const long double pz128 = 0.0L, nz128 = -0.0L;

   CHECK( std::signbit(Cl6(nz32)));
   CHECK(!std::signbit(Cl6(pz32)));
   CHECK( std::signbit(poly_Cl6(nz32)));
   CHECK(!std::signbit(poly_Cl6(pz32)));

   CHECK( std::signbit(Cl6(nz64)));
   CHECK(!std::signbit(Cl6(pz64)));
   CHECK( std::signbit(poly_Cl6(nz64)));
   CHECK(!std::signbit(poly_Cl6(pz64)));
#ifdef ENABLE_FORTRAN
   CHECK( std::signbit(poly_Cl6_fortran(nz64)));
   CHECK(!std::signbit(poly_Cl6_fortran(pz64)));
#endif

   CHECK( std::signbit(Cl6(nz128)));
   CHECK(!std::signbit(Cl6(pz128)));
   CHECK( std::signbit(poly_Cl6(nz128)));
   CHECK(!std::signbit(poly_Cl6(pz128)));
}

TEST_CASE("test_real_fixed_values")
{
   const auto eps64  = std::pow(10.0 , -std::numeric_limits<double>::digits10);
   const auto eps128 = std::pow(10.0L, -std::numeric_limits<long double>::digits10);
   const std::string filename(std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "Cl6.txt");
   const auto fixed_values = polylogarithm::test::read_reals_from_file<long double>(filename);

   for (auto v: fixed_values) {
      const auto x128 = v.first;
      const auto x64 = static_cast<double>(x128);
      const auto cl128_expected = v.second;
      const auto cl64_expected = static_cast<double>(cl128_expected);

      const auto cl64_li      = Cl6_via_Li6(x64);
      const auto cl64_poly    = polylogarithm::Cl6(x64);
      const auto cl64_poly_c  = cl6(x64);
#ifdef ENABLE_FORTRAN
      const auto cl64_poly_f  = poly_Cl6_fortran(x64);
#endif
      const auto cl128_li     = Cl6_via_Li6(x128);
      const auto cl128_poly   = polylogarithm::Cl6(x128);
      const auto cl128_poly_c = cl6l(x128);

      INFO("x(64)         = " << x64);
      INFO("Cl6(64)  real = " << cl64_expected  << " (expected)");
      INFO("Cl6(64)  real = " << cl64_poly      << " (polylogarithm C++)");
      INFO("Cl6(64)  real = " << cl64_poly_c    << " (polylogarithm C)");
#ifdef ENABLE_FORTRAN
      INFO("Cl6(64)  real = " << cl64_poly_f    << " (polylogarithm Fortran)");
#endif
      INFO("Cl6(64)  real = " << cl64_li        << " (via Li6 C++)");
      INFO("------------------------------------------------------------");
      INFO("x(128)        = " << x128);
      INFO("Cl6(128) real = " << cl128_expected << " (expected)");
      INFO("Cl6(128) real = " << cl128_poly     << " (polylogarithm C++)");
      INFO("Cl6(128) real = " << cl128_poly_c   << " (polylogarithm C)");
      INFO("Cl6(128) real = " << cl128_li       << " (via Li6 C++)");

      CHECK_CLOSE(cl64_li     , cl64_expected , 2*eps64 );
      CHECK_CLOSE(cl64_poly   , cl64_expected , 2*eps64 );
      CHECK_CLOSE(cl64_poly_c , cl64_expected , 2*eps64 );
#ifdef ENABLE_FORTRAN
      CHECK_CLOSE(cl64_poly_f , cl64_expected , 2*eps64 );
#endif
      CHECK_CLOSE(cl128_li    , cl128_expected, 4*eps128);
      CHECK_CLOSE(cl128_poly  , cl128_expected, 8*eps128);
      CHECK_CLOSE(cl128_poly_c, cl128_expected, 8*eps128);
   }
}
