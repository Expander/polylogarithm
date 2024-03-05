#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "Sl.hpp"
#include "Li.hpp"
#include "read_data.hpp"
#include "test.hpp"
#include <cmath>
#include <complex>

double Sl_via_Li(int64_t n, double x)
{
   const std::complex<double> li = polylogarithm::Li(n, std::polar(1.0, x));

   if (n % 2 == 0) {
      return std::real(li);
   }

   return std::imag(li);
}

// tests signbit for 0.0 and -0.0 arguments
TEST_CASE("test_signed_zero")
{
   // skip test if platform does not supprt signed zero
   if (!has_signed_zero()) {
      return;
   }

   using polylogarithm::Sl;

   const double pz64 = 0.0, nz64 = -0.0;

   for (int64_t n = 0; n <= 100; n++) {
      INFO("2*n + 1 = " << 2*n + 1);
      CHECK( std::signbit(Sl(2*n + 1, nz64)));
      CHECK(!std::signbit(Sl(2*n + 1, pz64)));
   }
}

TEST_CASE("test_fixed_values")
{
   const int ni[] = {
       1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
      11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
      21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
      31, 1000, 1001
   };

   for (const auto n: ni) {
      const std::string filename(std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "Sl" + std::to_string(n) + ".txt");
      const auto fixed_values = polylogarithm::test::read_reals_from_file<double>(filename);

      for (const auto& v: fixed_values) {
         const auto x = v.first;
         const auto cl_expected = v.second;
         INFO("n = " << n << ", x = " << x);
         CHECK_CLOSE(polylogarithm::Sl(n, x), cl_expected, 1e-13);
         CHECK_CLOSE(Sl_via_Li(n, x), cl_expected, 1e-9);
      }
   }
}
