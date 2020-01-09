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

   CHECK_CLOSE_COMPLEX(Li5(0), 0, 1e-15);
   CHECK_CLOSE_COMPLEX(Li5(1), zeta5, 1e-15);
   CHECK_CLOSE_COMPLEX(Li5(-1), -15.*zeta5/16., 1e-15);
   CHECK_CLOSE_COMPLEX(Li5(0.5), 0.5084005792422687, 1e-15);
}

TEST_CASE("test_fixed_values")
{
   const auto eps64  = 1e-14;  // @todo increase precision
   const std::string filename(std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "Li5.txt");
   const auto fixed_values = polylogarithm::test::read_from_file<long double>(filename);

   for (auto v: fixed_values) {
      const auto z128 = v.first;
      const auto z64 = to_c64(z128);
      const auto li128_expected = v.second;
      const auto li64_expected = to_c64(li128_expected);

      const auto li64_cmpl  = polylogarithm::Li5(z64);

      INFO("z(128)        = " << z128);
      INFO("Li3(64)  cmpl = " << li64_expected << " (expected)");
      INFO("Li3(128) cmpl = " << li128_expected << " (expected)");
      INFO("Li3(64)  cmpl = " << li64_cmpl << " (polylogarithm)");

      CHECK_CLOSE_COMPLEX(li64_cmpl , li64_expected , eps64);
   }
}
