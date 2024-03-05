#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "c_wrappers.h"
#include "Cl4.hpp"
#include "fortran_wrappers.h"
#include "Li4.hpp"
#include "read_data.hpp"
#include "test.hpp"
#include <cmath>
#include <complex>
#include <vector>

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

double poly_Cl4(double x) {
   return cl4(x);
}

long double poly_Cl4(long double x) {
   return cl4l(x);
}

#ifdef ENABLE_FORTRAN

double poly_Cl4_fortran(double x) {
   double res{};
   cl4_fortran(&x, &res);
   return res;
}

#endif

double Cl4_via_Li4(double x) noexcept
{
   return std::imag(polylogarithm::Li4(std::polar(1.0, x)));
}

long double Cl4_via_Li4(long double x) noexcept
{
   return std::imag(polylogarithm::Li4(std::polar(1.0L, x)));
}

TEST_CASE("test_special_values")
{
   using polylogarithm::Cl4;
   const double pi  = M_PI;

   // Cl_4(Pi/2) = DirichletBeta(4)
   CHECK_CLOSE(Cl4(pi/2.), 0.988944551741105336108422633228377821315860887062733910, 1e-15);
}

TEST_CASE("test_duplication_formula")
{
   using polylogarithm::Cl4;
   const double pi  = M_PI;

   const auto thetas = float_range(0., pi, 100);

   for (const auto t: thetas) {
      const auto rel = Cl4(2*t) - (8*Cl4(t) - 8*Cl4(pi - t));
      // converges not so good around t = 0, t = pi
      CHECK_SMALL(rel, 1e-8);
   }
}

TEST_CASE("test_roots")
{
   using polylogarithm::Cl4;
   const double pi  = M_PI;

   for (int k = -10; k < 10; k++) {
      CHECK_SMALL(Cl4(k*pi), 1e-10);
   }
}

// tests signbit for 0.0 and -0.0 arguments
TEST_CASE("test_signed_zero")
{
   // skip test if platform does not supprt signed zero
   if (!has_signed_zero()) {
      return;
   }

   using polylogarithm::Cl4;

   const float  pz32 = 0.0f, nz32 = -0.0f;
   const double pz64 = 0.0, nz64 = -0.0;
   const long double pz128 = 0.0L, nz128 = -0.0L;

   CHECK( std::signbit(Cl4(nz32)));
   CHECK(!std::signbit(Cl4(pz32)));
   CHECK( std::signbit(poly_Cl4(nz32)));
   CHECK(!std::signbit(poly_Cl4(pz32)));

   CHECK( std::signbit(Cl4(nz64)));
   CHECK(!std::signbit(Cl4(pz64)));
   CHECK( std::signbit(poly_Cl4(nz64)));
   CHECK(!std::signbit(poly_Cl4(pz64)));
#ifdef ENABLE_FORTRAN
   CHECK( std::signbit(poly_Cl4_fortran(nz64)));
   CHECK(!std::signbit(poly_Cl4_fortran(pz64)));
#endif

   CHECK( std::signbit(Cl4(nz128)));
   CHECK(!std::signbit(Cl4(pz128)));
   CHECK( std::signbit(poly_Cl4(nz128)));
   CHECK(!std::signbit(poly_Cl4(pz128)));
}

TEST_CASE("test_real_fixed_values")
{
   const auto eps64  = std::pow(10.0 , -std::numeric_limits<double>::digits10);
   const auto eps128 = std::pow(10.0L, -std::numeric_limits<long double>::digits10);
   const std::string filename(std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "Cl4.txt");
   const auto fixed_values = polylogarithm::test::read_reals_from_file<long double>(filename);

   for (auto v: fixed_values) {
      const auto x128 = v.first;
      const auto x64 = static_cast<double>(x128);
      const auto cl128_expected = v.second;
      const auto cl64_expected = static_cast<double>(cl128_expected);

      const auto cl64_li4     = Cl4_via_Li4(x64);
      const auto cl64_poly    = polylogarithm::Cl4(x64);
      const auto cl64_poly_c  = cl4(x64);
#ifdef ENABLE_FORTRAN
      const auto cl64_poly_f  = poly_Cl4_fortran(x64);
#endif
      const auto cl128_poly   = polylogarithm::Cl4(x128);
      const auto cl128_poly_c = cl4l(x128);
      const auto cl128_li4    = Cl4_via_Li4(x128);

      INFO("x(64)         = " << x64);
      INFO("Cl4(64)  real = " << cl64_expected  << " (expected)");
      INFO("Cl4(64)  real = " << cl64_poly      << " (polylogarithm C++)");
      INFO("Cl4(64)  real = " << cl64_poly_c    << " (polylogarithm C)");
#ifdef ENABLE_FORTRAN
      INFO("Cl4(64)  real = " << cl64_poly_f    << " (polylogarithm Fortran)");
#endif
      INFO("Cl4(64)  real = " << cl64_li4       << " (via Li4 C++)");
      INFO("------------------------------------------------------------");
      INFO("x(128)        = " << x128);
      INFO("Cl4(128) real = " << cl128_expected << " (expected)");
      INFO("Cl4(128) real = " << cl128_poly     << " (polylogarithm C++)");
      INFO("Cl4(128) real = " << cl128_poly_c   << " (polylogarithm C)");
      INFO("Cl4(128) real = " << cl128_li4      << " (via Li4 C++)");

      CHECK_CLOSE(cl64_li4    , cl64_expected , 2*eps64 );
      CHECK_CLOSE(cl64_poly   , cl64_expected , 2*eps64 );
      CHECK_CLOSE(cl64_poly_c , cl64_expected , 2*eps64 );
#ifdef ENABLE_FORTRAN
      CHECK_CLOSE(cl64_poly_f , cl64_expected , 2*eps64 );
#endif
      CHECK_CLOSE(cl128_li4   , cl128_expected, 2*eps128);
      CHECK_CLOSE(cl128_poly  , cl128_expected, 2*eps128);
      CHECK_CLOSE(cl128_poly_c, cl128_expected, 2*eps128);
   }
}
