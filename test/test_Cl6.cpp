#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "Li6.hpp"
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
   using namespace polylogarithm;
   const double pi  = M_PI;

   // Cl_6(Pi/2) = DirichletBeta(6)
   CHECK_CLOSE(Cl6(pi/2.), 0.998685222218438135441600787860206549678364546126514411, 1e-15);
}

TEST_CASE("test_duplication_formula")
{
   using namespace polylogarithm;
   const double pi  = M_PI;

   const auto thetas = float_range(0., pi, 100);

   for (const auto t: thetas) {
      const auto rel = Cl6(2*t) - (32*Cl6(t) - 32*Cl6(pi - t));
      CHECK_SMALL(rel, 1e-7);
   }
}

TEST_CASE("test_roots")
{
   using namespace polylogarithm;
   const double pi  = M_PI;

   for (int k = -10; k < 10; k++) {
      CHECK_SMALL(Cl6(k*pi), 1e-10);
   }
}
