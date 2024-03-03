#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "alt.h"
#include "c_wrappers.h"
#include "fortran_wrappers.h"
#include "Li6.hpp"
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

std::complex<double> poly_Li6(std::complex<double> z) {
   double re{}, im{};
   cli6_c(std::real(z), std::imag(z), &re, &im);
   return { re, im };
}

#ifdef ENABLE_FORTRAN

std::complex<double> poly_Li6_fortran(std::complex<double> z) {
   const double re = std::real(z);
   const double im = std::imag(z);
   double res_re{}, res_im{};
   cli6_fortran(&re, &im, &res_re, &res_im);
   return { res_re, res_im };
}

#endif

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

// tests signbit for 0.0 and -0.0 arguments
TEST_CASE("test_signed_zero")
{
   using polylogarithm::Li6;

   const float  pz32 = 0.0f, nz32 = -0.0f;
   const double pz64 = 0.0, nz64 = -0.0;
   const long double pz128 = 0.0L, nz128 = -0.0L;

   // complex Li6
   CHECK( std::signbit(std::real(Li6(std::complex<double>(nz64, nz64)))));
   CHECK( std::signbit(std::imag(Li6(std::complex<double>(nz64, nz64)))));
   CHECK(!std::signbit(std::real(Li6(std::complex<double>(pz64, nz64)))));
   CHECK( std::signbit(std::imag(Li6(std::complex<double>(pz64, nz64)))));
   CHECK( std::signbit(std::real(Li6(std::complex<double>(nz64, pz64)))));
   CHECK(!std::signbit(std::imag(Li6(std::complex<double>(nz64, pz64)))));
   CHECK(!std::signbit(std::real(Li6(std::complex<double>(pz64, pz64)))));
   CHECK(!std::signbit(std::imag(Li6(std::complex<double>(pz64, pz64)))));

   CHECK( std::signbit(std::real(Li6(std::complex<long double>(nz128, nz128)))));
   CHECK( std::signbit(std::imag(Li6(std::complex<long double>(nz128, nz128)))));
   CHECK(!std::signbit(std::real(Li6(std::complex<long double>(pz128, nz128)))));
   CHECK( std::signbit(std::imag(Li6(std::complex<long double>(pz128, nz128)))));
   CHECK( std::signbit(std::real(Li6(std::complex<long double>(nz128, pz128)))));
   CHECK(!std::signbit(std::imag(Li6(std::complex<long double>(nz128, pz128)))));
   CHECK(!std::signbit(std::real(Li6(std::complex<long double>(pz128, pz128)))));
   CHECK(!std::signbit(std::imag(Li6(std::complex<long double>(pz128, pz128)))));
}

template<typename T>
struct Data {
   std::complex<T> z;
   std::complex<T> li_expected;
   T eps;
};

TEST_CASE("test_overflow")
{
   using polylogarithm::Li6;

   const double eps64 = std::pow(10.0, -std::numeric_limits<double>::digits10);
   const long double eps128 = std::pow(10.0L, -std::numeric_limits<long double>::digits10);

   const Data<double> data64[] = {
      {{1e300, 1.0}, {-1.5086876165613597e14, 4.11768711823317e12}, 2*eps64},
      {{1.0, 1e300}, {-1.5090387516918862e14, 2058950021167.7976}, eps64}
   };

   for (const auto& d : data64) {
      CHECK_CLOSE_COMPLEX(Li6(d.z), d.li_expected, d.eps);
      CHECK_CLOSE_COMPLEX(poly_Li6(d.z), d.li_expected, d.eps);
#ifdef ENABLE_FORTRAN
      CHECK_CLOSE_COMPLEX(poly_Li6_fortran(d.z), d.li_expected, d.eps);
#endif
   }

   const Data<long double> data128[] = {
      {{1e4000L, 1.0L}, {-8.478540098238822566033321355282680660651e20L, 1.7351899625805641348387419061716509476e18L}, eps128}
   };

   for (const auto& d : data128) {
      CHECK_CLOSE_COMPLEX(Li6(d.z), d.li_expected, d.eps);
      CHECK_CLOSE_COMPLEX(poly_Li6(d.z), d.li_expected, d.eps);
   }
}

TEST_CASE("test_fixed_values")
{
   const auto eps64  = std::pow(10.0 , -std::numeric_limits<double>::digits10);
   const auto eps128 = std::pow(10.0L, -std::numeric_limits<long double>::digits10);

   const std::string filename(std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "Li6.txt");
   const auto fixed_values = polylogarithm::test::read_from_file<long double>(filename);

   for (auto v: fixed_values) {
      const auto z128 = v.first;
      const auto z64 = to<double>(z128);
      const auto li128_expected = v.second;
      const auto li64_expected = to<double>(li128_expected);

      const auto li64_cmpl  = polylogarithm::Li6(z64);
      const auto li128_cmpl = polylogarithm::Li6(z128);
      const auto li64_cmpl_c  = poly_Li6(z64);
#ifdef ENABLE_FORTRAN
      const auto li64_cmpl_f  = poly_Li6_fortran(z64);
#endif
      const auto li128_cmpl_c = poly_Li6(z128);

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
      CHECK_CLOSE_COMPLEX(li64_cmpl_f , li64_expected , 3*eps64);
#endif
      CHECK_CLOSE_COMPLEX(li128_cmpl  , li128_expected, 2*eps128);
      CHECK_CLOSE_COMPLEX(li128_cmpl_c, li128_expected, 2*eps128);
   }
}
