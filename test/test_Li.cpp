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
   using polylogarithm::Li;

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
         const auto li = Li(ni.n, z);
         INFO("n = " << ni.n << ", z = " << z);
         CHECK_CLOSE_COMPLEX(li, li_expected, ni.eps);
      }
   }

   // value close to boundary between series 1 and 2 in arXiv:2010.09860
   CHECK_CLOSE_COMPLEX(Li(-2, -0.50001), -0.074072592582716422, 1e-14);

   const auto z = std::complex<double>(1.5, 0.0);
   CHECK_CLOSE_COMPLEX(Li(10, z), std::complex<double>(1.5022603281703005298, -2.56429642116111388671e-9), 1e-14);
   CHECK_CLOSE_COMPLEX(Li(10, -z), std::complex<double>(-1.4978556954869267594, 0.0), 1e-14);
}

// tests signbit for 0.0 and -0.0 arguments
TEST_CASE("test_signed_zero")
{
   using polylogarithm::Li;

   const float  pz32 = 0.0f, nz32 = -0.0f;
   const double pz64 = 0.0, nz64 = -0.0;
   const long double pz128 = 0.0L, nz128 = -0.0L;

   // complex Li
   for (int n = -20; n <= 20; ++n) {
      CHECK( std::signbit(std::real(Li(n, std::complex<double>(nz64, nz64)))));
      CHECK( std::signbit(std::imag(Li(n, std::complex<double>(nz64, nz64)))));
      CHECK(!std::signbit(std::real(Li(n, std::complex<double>(pz64, nz64)))));
      CHECK( std::signbit(std::imag(Li(n, std::complex<double>(pz64, nz64)))));
      CHECK( std::signbit(std::real(Li(n, std::complex<double>(nz64, pz64)))));
      CHECK(!std::signbit(std::imag(Li(n, std::complex<double>(nz64, pz64)))));
      CHECK(!std::signbit(std::real(Li(n, std::complex<double>(pz64, pz64)))));
      CHECK(!std::signbit(std::imag(Li(n, std::complex<double>(pz64, pz64)))));
   }
}

template<typename T>
struct Data {
   int n;
   std::complex<T> z;
   std::complex<T> li_expected;
   T eps;
};

TEST_CASE("test_overflow")
{
   using polylogarithm::Li;

   const double eps64 = std::pow(10.0, -std::numeric_limits<double>::digits10);

   const Data<double> data64[] = {
      {7, {1e300, 1.0}, {-1.4886831990993457e16, 4.74066248802866e14}, eps64},
      {7, {1.0, 1e300}, {-1.489168315226607e16, 2.3705150998401e14}, eps64}
   };

   for (const auto& d : data64) {
      CHECK_CLOSE_COMPLEX(Li(d.n, d.z), d.li_expected, d.eps);
   }
}
