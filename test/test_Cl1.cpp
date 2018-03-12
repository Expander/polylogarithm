#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "Li1.hpp"
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

TEST_CASE("test_duplication_formula")
{
   using namespace polylogarithm;
   const double pi  = M_PI;

   const auto thetas = float_range(0. + 0.01, pi, 100);

   for (const auto t: thetas) {
      const auto rel = Cl1(2*t) - (Cl1(t) + Cl1(pi - t));
      CHECK_SMALL(rel, 1e-10);
   }
}
