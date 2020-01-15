#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "Li5.hpp"
#include "read_data.hpp"
#include <cmath>
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

TEST_CASE("test_special_values")
{
   using polylogarithm::Li5;

   const double zeta5 = 1.0369277551433699;
   const std::complex<double> zero(0.0, 0.0);
   const std::complex<double> one(1.0, 0.0);
   const std::complex<double> mone(-1.0, 0.0);
   const std::complex<double> half(0.5, 0.0);

   CHECK_CLOSE_COMPLEX(Li5(zero), 0, 1e-15);
   CHECK_CLOSE_COMPLEX(Li5(one), zeta5, 1e-15);
   CHECK_CLOSE_COMPLEX(Li5(mone), -15.*zeta5/16., 1e-15);
   CHECK_CLOSE_COMPLEX(Li5(half), 0.5084005792422687, 1e-15);
}

TEST_CASE("test_fixed_values")
{
   const std::string filename(std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "Li5.txt");
   const auto fixed_values = polylogarithm::test::read_from_file<long double>(filename);

   for (auto v: fixed_values) {
      const auto z128 = v.first;
      const auto z64 = to_c64(z128);
      const auto li128_expected = v.second;
      const auto li64_expected = to_c64(li128_expected);

      const auto li64_cmpl  = polylogarithm::Li5(z64);
      const auto li128_cmpl = polylogarithm::Li5(z128);

      INFO("z(128)        = " << z128);
      INFO("Li5(64)  cmpl = " << li64_expected << " (expected)");
      INFO("Li5(128) cmpl = " << li128_expected << " (expected)");
      INFO("Li5(64)  cmpl = " << li64_cmpl << " (polylogarithm)");
      INFO("Li5(128) cmpl = " << li128_cmpl << " (polylogarithm)");

      CHECK_CLOSE_COMPLEX(li64_cmpl , li64_expected , 3e-15);
      CHECK_CLOSE_COMPLEX(li128_cmpl, li128_expected, 2e-18);
   }
}
