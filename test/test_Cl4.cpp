#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "Cl4.hpp"
#include "Li4.hpp"
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
      const auto cl128_poly   = polylogarithm::Cl4(x128);
      const auto cl128_li4    = Cl4_via_Li4(x128);

      INFO("x(64)         = " << x64);
      INFO("Cl4(64)  real = " << cl64_expected  << " (expected)");
      INFO("Cl4(64)  real = " << cl64_poly      << " (polylogarithm C++)");
      INFO("Cl4(64)  real = " << cl64_li4       << " (via Li4 C++)");
      INFO("------------------------------------------------------------");
      INFO("x(128)        = " << x128);
      INFO("Cl4(128) real = " << cl128_expected << " (expected)");
      INFO("Cl4(128) real = " << cl128_poly     << " (polylogarithm C++)");
      INFO("Cl4(128) real = " << cl128_li4      << " (via Li4 C++)");

      CHECK_CLOSE(cl64_poly , cl64_expected , 2*eps64 );
      CHECK_CLOSE(cl128_poly, cl128_expected, 2*eps128);
   }
}
