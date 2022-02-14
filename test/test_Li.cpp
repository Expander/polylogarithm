#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "Li.hpp"
#include "read_data.hpp"
#include <cmath>
#include <complex>
#include <string>

#define CHECK_CLOSE(a,b,eps) CHECK((a) == doctest::Approx(b).epsilon(eps))

#define CHECK_CLOSE_COMPLEX(a,b,eps) do {                               \
      CHECK_CLOSE(std::real(a), std::real(b), (eps));                   \
      CHECK_CLOSE(std::imag(a), std::imag(b), (eps));                   \
   } while (0)

TEST_CASE("test_complex_fixed_values")
{
   const auto eps = 1e-13;
   const int ni[] = {-2, -1, 0, 1, 2, 3, 4, 5, 6, 100};

   for (const auto n: ni) {
      const std::string filename(std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "Li" + std::to_string(n) + ".txt");
      const auto values = polylogarithm::test::read_from_file<double>(filename);

      for (auto v: values) {
         const auto z = v.first;
         const auto li_expected = v.second;
         INFO("n = " << n << ", z = " << z);
         CHECK_CLOSE_COMPLEX(polylogarithm::Li(n,z), li_expected, eps);
      }
   }
}
