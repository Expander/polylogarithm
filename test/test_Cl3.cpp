#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "alt.h"
#include "c_wrappers.h"
#include "Cl3.hpp"
#include "fortran_wrappers.h"
#include "Li3.hpp"
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

#ifdef ENABLE_FORTRAN

double poly_Cl3_fortran(double x) {
   double res{};
   cl3_fortran(&x, &res);
   return res;
}

#endif

double Cl3_via_Li3(double x) noexcept
{
   return std::real(polylogarithm::Li3(std::polar(1.0, x)));
}

long double Cl3_via_Li3(long double x) noexcept
{
   return std::real(polylogarithm::Li3(std::polar(1.0L, x)));
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

      const auto cl64_li3     = Cl3_via_Li3(x64);
      const auto cl64_pade    = clausen_3_pade(x64);
      const auto cl64_poly    = polylogarithm::Cl3(x64);
      const auto cl64_poly_c  = cl3(x64);
#ifdef ENABLE_FORTRAN
      const auto cl64_poly_f  = poly_Cl3_fortran(x64);
#endif
      const auto cl64_wu      = clausen_3_wu(x64);
      const auto cl128_poly   = polylogarithm::Cl3(x128);
      const auto cl128_poly_c = cl3l(x128);
      const auto cl128_li3    = Cl3_via_Li3(x128);

      INFO("x(64)         = " << x64);
      INFO("Cl3(64)  real = " << cl64_expected  << " (expected)");
      INFO("Cl3(64)  real = " << cl64_poly      << " (polylogarithm C++)");
      INFO("Cl3(64)  real = " << cl64_poly_c    << " (polylogarithm C)");
#ifdef ENABLE_FORTRAN
      INFO("Cl3(64)  real = " << cl64_poly_f    << " (polylogarithm Fortran)");
#endif
      INFO("Cl3(64)  real = " << cl64_li3       << " (via Li3 C++)");
      INFO("Cl3(64)  real = " << cl64_pade      << " (Pade)");
      INFO("Cl3(64)  real = " << cl64_wu        << " (Wu et.al.)");
      INFO("------------------------------------------------------------");
      INFO("x(128)        = " << x128);
      INFO("Cl3(128) real = " << cl128_expected << " (expected)");
      INFO("Cl3(128) real = " << cl128_poly     << " (polylogarithm C++)");
      INFO("Cl3(128) real = " << cl128_poly_c   << " (polylogarithm C)");
      INFO("Cl3(128) real = " << cl128_li3      << " (via Li3 C++)");

      CHECK_CLOSE(cl64_li3    , cl64_expected , 2*eps64 );
      CHECK_CLOSE(cl64_poly   , cl64_expected , 2*eps64 );
      CHECK_CLOSE(cl64_poly_c , cl64_expected , 2*eps64 );
#ifdef ENABLE_FORTRAN
      CHECK_CLOSE(cl64_poly_f , cl64_expected , 2*eps64 );
#endif
      CHECK_CLOSE(cl64_pade   , cl64_expected , 2*eps64 );
      CHECK_CLOSE(cl64_wu     , cl64_expected , 2*eps64 );
      CHECK_CLOSE(cl128_poly  , cl128_expected, 5*eps128);
      CHECK_CLOSE(cl128_poly_c, cl128_expected, 5*eps128);
      CHECK_CLOSE(cl128_li3   , cl128_expected, 2*eps128);

      // test symmetries
      if (std::abs(std::fmod(x64, 2*M_PI)) > 0.1 && std::abs(x64 - 2*M_PI) > 0.1) {
         CHECK_CLOSE(  polylogarithm::Cl3(x64 + 2*M_PI),  cl64_expected ,  10*eps64);
         CHECK_CLOSE(  polylogarithm::Cl3(x64 - 2*M_PI),  cl64_expected ,  10*eps64);
         CHECK_CLOSE(  polylogarithm::Cl3(-x64        ),  cl64_expected ,  10*eps64);
         CHECK_CLOSE(  polylogarithm::Cl3(-x64        ),  cl64_expected ,  10*eps64);

         CHECK_CLOSE(                 cl3(x64 + 2*M_PI),  cl64_expected ,  10*eps64);
         CHECK_CLOSE(                 cl3(x64 - 2*M_PI),  cl64_expected ,  10*eps64);
         CHECK_CLOSE(                 cl3(-x64        ),  cl64_expected ,  10*eps64);
         CHECK_CLOSE(                 cl3(-x64        ),  cl64_expected ,  10*eps64);

#ifdef ENABLE_FORTRAN
         CHECK_CLOSE(    poly_Cl3_fortran(x64 + 2*M_PI),  cl64_expected ,  10*eps64);
         CHECK_CLOSE(    poly_Cl3_fortran(x64 - 2*M_PI),  cl64_expected ,  10*eps64);
         CHECK_CLOSE(    poly_Cl3_fortran(-x64        ),  cl64_expected ,  10*eps64);
         CHECK_CLOSE(    poly_Cl3_fortran(-x64        ),  cl64_expected ,  10*eps64);
#endif

         CHECK_CLOSE(polylogarithm::Cl3(x128 + 2*M_PIL),  cl128_expected, 10*eps128);
         CHECK_CLOSE(polylogarithm::Cl3(x128 - 2*M_PIL),  cl128_expected, 10*eps128);
         CHECK_CLOSE(polylogarithm::Cl3(-x128         ),  cl128_expected, 10*eps128);
         CHECK_CLOSE(polylogarithm::Cl3(-x128         ),  cl128_expected, 10*eps128);
      }
   }
}
