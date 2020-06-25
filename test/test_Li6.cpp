#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "alt.h"
#include "Li6.hpp"
#include "polylogarithm.h"
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

std::complex<double> to_c64(std::complex<long double> z)
{
   return std::complex<double>(std::real(z), std::imag(z));
}

std::complex<double> poly_Li6(std::complex<double> z) {
   double re{}, im{};
   cli6_c(std::real(z), std::imag(z), &re, &im);
   return { re, im };
}

std::complex<long double> poly_Li6(std::complex<long double> z) {
   long double re{}, im{};
   cli6l_c(std::real(z), std::imag(z), &re, &im);
   return { re, im };
}

TEST_CASE("test_special_values")
{
   using polylogarithm::Li6;

   const double eps = std::pow(10.0, -std::numeric_limits<double>::digits10);
   const double zeta6 = 1.0173430619844491397145179297909;
   const std::complex<double> zero(0.0, 0.0);
   const std::complex<double> one(1.0, 0.0);
   const std::complex<double> mone(-1.0, 0.0);
   const std::complex<double> half(0.5, 0.0);

   CHECK_CLOSE_COMPLEX(Li6(zero), 0, eps);
   CHECK_CLOSE_COMPLEX(Li6(one), zeta6, eps);
   CHECK_CLOSE_COMPLEX(Li6(mone), -31.*zeta6/32., eps);
   CHECK_CLOSE_COMPLEX(Li6(half), 0.5040953978039886, eps);
}

TEST_CASE("test_fixed_values")
{
   const auto eps64  = std::pow(10.0 , -std::numeric_limits<double>::digits10);
   const auto eps128 = std::pow(10.0L, -std::numeric_limits<long double>::digits10);

   const std::string filename(std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "Li6.txt");
   const auto fixed_values = polylogarithm::test::read_from_file<long double>(filename);

   for (auto v: fixed_values) {
      const auto z128 = v.first;
      const auto z64 = to_c64(z128);
      const auto li128_expected = v.second;
      const auto li64_expected = to_c64(li128_expected);

      const auto li64_cmpl  = polylogarithm::Li6(z64);
      const auto li128_cmpl = polylogarithm::Li6(z128);
      const auto li64_cmpl_c  = poly_Li6(z64);
      const auto li128_cmpl_c = poly_Li6(z128);

      INFO("z(128)        = " << z128);
      INFO("Li4(64)  cmpl = " << li64_expected  << " (expected)");
      INFO("Li4(64)  cmpl = " << li64_cmpl      << " (polylogarithm C++)");
      INFO("Li4(64)  cmpl = " << li64_cmpl_c    << " (polylogarithm C)");
      INFO("Li4(128) cmpl = " << li128_expected << " (expected)");
      INFO("Li4(128) cmpl = " << li128_cmpl     << " (polylogarithm C++)");
      INFO("Li4(128) cmpl = " << li128_cmpl_c   << " (polylogarithm C)");

      CHECK_CLOSE_COMPLEX(li64_cmpl   , li64_expected , 2*eps64);
      CHECK_CLOSE_COMPLEX(li64_cmpl_c , li64_expected , 2*eps64);
      CHECK_CLOSE_COMPLEX(li128_cmpl  , li128_expected, 2*eps128);
      CHECK_CLOSE_COMPLEX(li128_cmpl_c, li128_expected, 2*eps128);
   }
}
