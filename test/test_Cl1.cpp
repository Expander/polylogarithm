#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "Cl1.hpp"
#include "read_data.hpp"
#include <cmath>

#define CHECK_CLOSE(a,b,eps) CHECK((a) == doctest::Approx(b).epsilon(eps))

TEST_CASE("test_real_fixed_values")
{
   const auto eps64  = std::pow(10.0 , -std::numeric_limits<double>::digits10);
   const std::string filename(std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "Cl1.txt");
   const auto fixed_values = polylogarithm::test::read_reals_from_file<long double>(filename);

   for (auto v: fixed_values) {
      const auto x128 = v.first;
      const auto x64 = static_cast<double>(x128);
      const auto cl128_expected = v.second;
      const auto cl64_expected = static_cast<double>(cl128_expected);

      const auto cl64_poly    = polylogarithm::Cl1(x64);

      INFO("x(64)         = " << x64);
      INFO("Cl1(64)  real = " << cl64_expected  << " (expected)");
      INFO("Cl1(64)  real = " << cl64_poly      << " (polylogarithm C++)");

      CHECK_CLOSE(cl64_poly, cl64_expected, 4*eps64);
   }
}
