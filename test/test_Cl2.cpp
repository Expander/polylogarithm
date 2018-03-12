#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "Li2.hpp"
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
   const double catalan = 0.91596559417721901505460351493238411077414937428167;

   CHECK_CLOSE(Cl2(pi/2.), catalan, 1e-15);
}

TEST_CASE("test_kummer_relation")
{
   using namespace polylogarithm;
   using std::exp;
   const double pi  = M_PI;
   const double z2  = 1.644934066848226436472415166646025189218949901206798437735558229;
   const std::complex<double> i(0.,1.);

   const auto thetas = float_range(0., 2*pi, 100);

   for (const auto t: thetas) {
      const auto lhs = Li2(exp(i*t));
      const auto rhs = z2 - t*(2*pi - t)/4. + i*Cl2(t);

      CHECK_CLOSE(std::real(lhs), std::real(rhs), 1e-15);
      CHECK_CLOSE(std::imag(lhs), std::imag(rhs), 1e-15);
   }
}

TEST_CASE("test_roots")
{
   using namespace polylogarithm;
   const double pi  = M_PI;

   for (int k = -10; k < 10; k++) {
      CHECK_SMALL(Cl2(k*pi), 1e-10);
   }
}
