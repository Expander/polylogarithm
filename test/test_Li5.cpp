#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "alt.h"
#include "c_wrappers.h"
#include "fortran_wrappers.h"
#include "Li5.hpp"
#include "read_data.hpp"
#include <cmath>
#include <limits>
#include <utility>

#define CHECK_CLOSE(a,b,eps) CHECK((a) == doctest::Approx(b).epsilon(eps))
#define CHECK_CLOSE_COMPLEX(a,b,eps) do {                               \
      CHECK_CLOSE(std::real(a), std::real(b), (eps));                   \
      CHECK_CLOSE(std::imag(a), std::imag(b), (eps));                   \
   } while (0)
#define CHECK_SMALL(a,eps) CHECK(std::abs(a) < (eps))

template <typename T, typename U>
std::complex<T> to(std::complex<U> z)
{
   return std::complex<T>(static_cast<T>(std::real(z)), static_cast<T>(std::imag(z)));
}

std::complex<double> poly_Li5(std::complex<double> z) {
   double re{}, im{};
   cli5_c(std::real(z), std::imag(z), &re, &im);
   return { re, im };
}

#ifdef ENABLE_FORTRAN

std::complex<double> poly_Li5_fortran(std::complex<double> z) {
   const double re = std::real(z);
   const double im = std::imag(z);
   double res_re{}, res_im{};
   cli5_fortran(&re, &im, &res_re, &res_im);
   return { res_re, res_im };
}

#endif

std::complex<long double> poly_Li5(std::complex<long double> z) {
   long double re{}, im{};
   cli5l_c(std::real(z), std::imag(z), &re, &im);
   return { re, im };
}

TEST_CASE("test_special_values")
{
   using polylogarithm::Li5;

   const double eps = std::pow(10.0, -std::numeric_limits<double>::digits10);
   const double zeta5 = 1.0369277551433699;
   const std::complex<double> zero(0.0, 0.0);
   const std::complex<double> one(1.0, 0.0);
   const std::complex<double> mone(-1.0, 0.0);
   const std::complex<double> half(0.5, 0.0);

   CHECK_CLOSE_COMPLEX(Li5(zero), 0, eps);
   CHECK_CLOSE_COMPLEX(Li5(one), zeta5, eps);
   CHECK_CLOSE_COMPLEX(Li5(mone), -15.*zeta5/16., eps);
   CHECK_CLOSE_COMPLEX(Li5(half), 0.5084005792422687, eps);

   {
      // value that cause overflow when squared
      const std::complex<double> z(1e300, 1.0);
      const std::complex<double> ze(-1.3105197831948743e12, -2.980481322754618e10);
      CHECK_CLOSE(std::real(Li5(z)), std::real(ze), eps);
      CHECK_CLOSE(std::real(poly_Li5(z)), std::real(ze), eps);
#ifdef ENABLE_FORTRAN
      CHECK_CLOSE(std::real(poly_Li5_fortran(z)), std::real(ze), eps);
#endif
   }

   {
      // values that cause overflow when squared
      const std::complex<double> z(1.0, 1e300);
      const std::complex<double> ze(-1.31072310968392418e12, 1.490286896860219e10);
      CHECK_CLOSE(std::real(Li5(z)), std::real(ze), eps);
      CHECK_CLOSE(std::real(poly_Li5(z)), std::real(ze), eps);
#ifdef ENABLE_FORTRAN
      CHECK_CLOSE(std::real(poly_Li5_fortran(z)), std::real(ze), eps);
#endif
   }

   {
      // value that cause overflow when squared
      const std::complex<long double> z(1e4000L, 1.0L);
      const std::complex<long double> ze(-5.523276910915025884712078063195501148879e17L, -9.419792822533112833122424859594537172e14L);
      const auto eps = std::pow(10.0L, -std::numeric_limits<long double>::digits10);
      CHECK_CLOSE(std::real(Li5(z)), std::real(ze), eps);
      CHECK_CLOSE(std::real(poly_Li5(z)), std::real(ze), eps);
   }
}

TEST_CASE("test_fixed_values")
{
   const auto eps64  = std::pow(10.0 , -std::numeric_limits<double>::digits10);
   const auto eps128 = std::pow(10.0L, -std::numeric_limits<long double>::digits10);

   const std::string filename(std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "Li5.txt");
   const auto fixed_values = polylogarithm::test::read_from_file<long double>(filename);

   for (auto v: fixed_values) {
      const auto z128 = v.first;
      const auto z64 = to<double>(z128);
      const auto li128_expected = v.second;
      const auto li64_expected = to<double>(li128_expected);

      const auto li64_cmpl  = polylogarithm::Li5(z64);
      const auto li128_cmpl = polylogarithm::Li5(z128);
      const auto li64_cmpl_c  = poly_Li5(z64);
#ifdef ENABLE_FORTRAN
      const auto li64_cmpl_f  = poly_Li5_fortran(z64);
#endif
      const auto li128_cmpl_c = poly_Li5(z128);

      INFO("z(128)        = " << z128);
      INFO("Li4(64)  cmpl = " << li64_expected  << " (expected)");
      INFO("Li4(64)  cmpl = " << li64_cmpl      << " (polylogarithm C++)");
      INFO("Li4(64)  cmpl = " << li64_cmpl_c    << " (polylogarithm C)");
#ifdef ENABLE_FORTRAN
      INFO("Li4(64)  cmpl = " << li64_cmpl_f    << " (polylogarithm Fortran)");
#endif
      INFO("Li4(128) cmpl = " << li128_expected << " (expected)");
      INFO("Li4(128) cmpl = " << li128_cmpl     << " (polylogarithm C++)");
      INFO("Li4(128) cmpl = " << li128_cmpl_c   << " (polylogarithm C)");

      CHECK_CLOSE_COMPLEX(li64_cmpl   , li64_expected , 2*eps64);
      CHECK_CLOSE_COMPLEX(li64_cmpl_c , li64_expected , 2*eps64);
#ifdef ENABLE_FORTRAN
      CHECK_CLOSE_COMPLEX(li64_cmpl_f , li64_expected , 2*eps64);
#endif
      CHECK_CLOSE_COMPLEX(li128_cmpl  , li128_expected, 2*eps128);
      CHECK_CLOSE_COMPLEX(li128_cmpl_c, li128_expected, 2*eps128);
   }
}
