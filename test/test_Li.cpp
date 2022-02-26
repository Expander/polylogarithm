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

/*

TEST_CASE("test_infinite_values")
{
   using polylogarithm::Li;

   for (int n = -10; n <= 0; ++n) {
      CHECK(std::isinf(std::real(Li(n, 1.0))));
      CHECK(std::isinf(std::imag(Li(n, 1.0))));
   }
}

*/

TEST_CASE("test_complex_fixed_values")
{
   const struct {
      int n;
      double eps;
   } nis[] = {
      // { -100, 1e-05 },
      {  -10, 1e-09 },
      {   -9, 1e-10 },
      {   -8, 1e-10 },
      {   -7, 1e-11 },
      {   -6, 1e-11 },
      {   -5, 1e-10 },
      {   -4, 1e-13 },
      {   -3, 1e-13 },
      {   -2, 1e-13 },
      {   -1, 1e-14 },
      {    0, 1e-14 },
      {    1, 1e-14 },
      {    2, 1e-14 },
      {    3, 1e-14 },
      {    4, 1e-14 },
      {    5, 1e-14 },
      {    6, 1e-14 },
      {  100, 1e-14 }
   };

   for (const auto ni: nis) {
      const std::string filename(std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "Li" + std::to_string(ni.n) + ".txt");
      const auto values = polylogarithm::test::read_from_file<double>(filename);

      for (auto v: values) {
         const auto z = v.first;
         const auto li_expected = v.second;
         const auto li = polylogarithm::Li(ni.n, z);
         INFO("n = " << ni.n << ", z = " << z);
         CHECK_CLOSE_COMPLEX(li, li_expected, ni.eps);
      }
   }

   // value close to boundary between series 1 and 2 in arXiv:2010.09860
   CHECK_CLOSE_COMPLEX(polylogarithm::Li(-2, -0.50001), -0.074072592582716422, 1e-14);
}
