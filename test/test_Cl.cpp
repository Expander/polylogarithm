#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "Li.hpp"
#include "Li1.hpp"
#include "Li2.hpp"
#include "Li3.hpp"
#include "Li4.hpp"
#include "Li5.hpp"
#include "Li6.hpp"
#include <cmath>
#include <complex>
#include <vector>

#define CHECK_CLOSE(a,b,eps) do {                       \
      if (std::isinf(a) && std::isinf(b))               \
         CHECK(true);                                   \
      else                                              \
         CHECK((a) == doctest::Approx(b).epsilon(eps)); \
   } while (0);
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

   CHECK_CLOSE(Cl(2,pi/2.), catalan, 1e-15);
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
      const auto lhs = Li(2,exp(i*t));
      const auto rhs = z2 - t*(2*pi - t)/4. + i*Cl(2,t);

      CHECK_CLOSE(std::real(lhs), std::real(rhs), 1e-15);
      CHECK_CLOSE(std::imag(lhs), std::imag(rhs), 1e-15);
   }
}

TEST_CASE("test_fixed_implementations")
{
   using namespace polylogarithm;

   const double pi = M_PI;
   const double eps = 1e-10;
   const auto thetas = float_range(0., 2*pi, 100);

   for (const auto t: thetas) {
      const auto cl1 = Cl1(t);
      const auto cl2 = Cl2(t);
      const auto cl3 = Cl3(t);
      const auto cl4 = Cl4(t);
      const auto cl5 = Cl5(t);
      const auto cl6 = Cl6(t);

      CHECK_CLOSE(cl1, Cl(1,t), eps);
      CHECK_CLOSE(cl2, Cl(2,t), eps);
      CHECK_CLOSE(cl3, Cl(3,t), eps);
      CHECK_CLOSE(cl4, Cl(4,t), eps);
      CHECK_CLOSE(cl5, Cl(5,t), eps);
      CHECK_CLOSE(cl6, Cl(6,t), eps);
   }
}

TEST_CASE("test_duplication_formula")
{
   using namespace polylogarithm;

   const double pi = M_PI;
   const double eps = 1e-9;
   const auto thetas = float_range(0., pi, 100);

   for (const auto t: thetas) {
      for (int m = 1; m < 20; m++) {
         const int sgn = m % 2 == 0 ? 1 : -1;
         const double lhs = Cl(m+1,2*t);
         const double rhs = std::pow(2,m)*(Cl(m+1,t) + sgn*Cl(m+1,pi-t));

         CHECK_CLOSE(lhs, rhs, eps);
      }
   }
}

TEST_CASE("test_roots")
{
   using namespace polylogarithm;
   const double pi  = M_PI;

   for (int k = -10; k < 10; k++) {
      CHECK_SMALL(Cl(2,k*pi), 1e-10);
   }
}
