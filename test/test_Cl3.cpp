#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "Li3.hpp"
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
   const double pi = M_PI;
   const double z3 = 1.202056903159594;

   // Cl_3(Pi/2) = -(3 Zeta[3])/32
   CHECK_CLOSE(Cl3(pi/2.), -3*z3/32., 1e-15);
   CHECK_CLOSE(Cl3(0.), z3, 1e-15);
   CHECK_CLOSE(Cl3(2*pi), z3, 1e-15);
}

TEST_CASE("test_duplication_formula")
{
   using namespace polylogarithm;
   const double pi  = M_PI;

   const auto thetas = float_range(0., pi, 100);

   for (const auto t: thetas) {
      const auto rel = Cl3(2*t) - (4*Cl3(t) + 4*Cl3(pi - t));
      CHECK_SMALL(rel, 2e-10);
   }
}
