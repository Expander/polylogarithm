#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "alt.h"
#include "Cl3.hpp"
#include "read_data.hpp"
#include <cmath>
#include <complex>
#include <vector>

#define CHECK_CLOSE(a,b,eps) CHECK((a) == doctest::Approx(b).epsilon(eps))
#define CHECK_SMALL(a,eps) CHECK(std::abs(a) <= (eps))

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
   using polylogarithm::Cl3;
   const double pi = M_PI;
   const double z3 = 1.202056903159594;

   // Cl_3(Pi/2) = -(3 Zeta[3])/32
   CHECK_CLOSE(Cl3(pi/2.), -3*z3/32., 1e-15);
   CHECK_CLOSE(Cl3(0.), z3, 1e-15);
   CHECK_CLOSE(Cl3(2*pi), z3, 1e-15);
}

TEST_CASE("test_duplication_formula")
{
   using polylogarithm::Cl3;
   const double pi = M_PI;

   const auto thetas = float_range(0., pi, 100);

   for (const auto t: thetas) {
      const auto rel = Cl3(2*t) - (4*Cl3(t) + 4*Cl3(pi - t));
      CHECK_SMALL(rel, 2e-10);
   }
}

TEST_CASE("test_real_fixed_values")
{
   const auto eps64  = std::pow(10.0 , -std::numeric_limits<double>::digits10);
   const auto eps128 = std::pow(10.0L, -std::numeric_limits<long double>::digits10);
   const std::string filename(std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "Cl3.txt");
   const auto fixed_values = polylogarithm::test::read_reals_from_file<long double>(filename);

   for (auto v: fixed_values) {
      const auto x128 = v.first;
      const auto x64 = static_cast<double>(x128);
      const auto cl128_expected = v.second;
      const auto cl64_expected = static_cast<double>(cl128_expected);

      const auto cl64_poly    = polylogarithm::Cl3(x64);
      const auto cl64_wu      = clausen_3_wu(x64);
      const auto cl128_poly   = polylogarithm::Cl3(x128);

      INFO("x(64)         = " << x64);
      INFO("Cl3(64)  real = " << cl64_expected  << " (expected)");
      INFO("Cl3(64)  real = " << cl64_poly      << " (polylogarithm C++)");
      INFO("Cl3(64)  real = " << cl64_wu        << " (Wu et.al.)");
      INFO("------------------------------------------------------------");
      INFO("x(128)        = " << x128);
      INFO("Cl3(128) real = " << cl128_expected << " (expected)");
      INFO("Cl3(128) real = " << cl128_poly     << " (polylogarithm C++)");

      CHECK_CLOSE(cl64_poly , cl64_expected , 2*eps64 );
      CHECK_CLOSE(cl64_wu   , cl64_expected , 2*eps64 );
      CHECK_CLOSE(cl128_poly, cl128_expected, 2*eps128);
   }
}
