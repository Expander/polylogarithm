#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "alt.h"
#include "c_wrappers.h"
#include "fortran_wrappers.h"
#include "Li3.hpp"
#include "bench.hpp"
#include "read_data.hpp"
#include "test.hpp"
#include <cmath>
#include <limits>
#include <utility>

std::complex<double> clog(std::complex<double> z) {
   std::complex<double> zf(z);
   // convert -0.0 to 0.0
   if (std::real(zf) == 0.0) { zf.real(0.0); }
   if (std::imag(zf) == 0.0) { zf.imag(0.0); }
   return std::log(zf);
}

double poly_Li3(double z) {
   return li3(z);
}

std::complex<double> poly_Li3(std::complex<double> z) {
   double re{}, im{};
   cli3_c(std::real(z), std::imag(z), &re, &im);
   return { re, im };
}

#ifdef ENABLE_FORTRAN

double poly_Li3_fortran(double x) {
   double res{};
   li3_fortran(&x, &res);
   return res;
}

std::complex<double> poly_Li3_fortran(std::complex<double> z) {
   const double re = std::real(z);
   const double im = std::imag(z);
   double res_re{}, res_im{};
   cli3_fortran(&re, &im, &res_re, &res_im);
   return { res_re, res_im };
}

#endif

std::complex<long double> poly_Li3(std::complex<long double> z) {
   long double re{}, im{};
   cli3l_c(std::real(z), std::imag(z), &re, &im);
   return { re, im };
}

std::complex<long double> tsil_Li3(std::complex<long double> z)
{
   long double re{}, im{};
   tsil_trilog_complex(std::real(z), std::imag(z), &re, &im);
   return { re, im };
}

const auto Relation_1 = [](std::complex<double> z) {
   using polylogarithm::Li3;
   return Li3(z) + Li3(-z) - Li3(z*z)/4.;
};

const auto Relation_2 = [](std::complex<double> z) -> std::complex<double> {
   using polylogarithm::Li3;
   using std::log;

   if (std::abs(z) == 0) {
      return {0.,0.};
   }

   // Relation does not seem to hold for Li[3,z], not even for
   // Mathematica's PolyLog[3,z], when 0 < Re[z] < 1.
   if (std::real(z) > 0. && std::real(z) < 1.) {
      return {0.,0.};
   }

   return Li3(z) - Li3(1./z) - (-pow3(clog(-z))/6. - M_PI*M_PI/6.*clog(-z));
};

const auto Relation_3 = [](std::complex<double> z) -> std::complex<double> {
   using polylogarithm::Li3;
   using std::log;
   const double zeta3 = 1.202056903159594;

   if (std::abs(std::real(1. - z)) < 1e-10) {
      return {0.,0.};
   }

   // Relation does not seem to hold for Li[3,z], not even for
   // Mathematica's PolyLog[3,z], when Re[z] < 0 and Im[z] = 0.
   if (std::real(z) <= 0. && std::imag(z) == 0.) {
      return {0.,0.};
   }

   return Li3(z) + Li3(1.-z) + Li3(1.-1./z)
      - (zeta3 + pow3(log(z))/6. + M_PI*M_PI/6.*clog(z) - 0.5*sqr(clog(z))*clog(1.-z));
};

TEST_CASE("test_special_values")
{
   using polylogarithm::Li3;

   const double ln2 = std::log(2.);
   const double pi2 = sqr(M_PI);
   const double zeta3 = 1.2020569031595942853997381615114;
   const double phi = 0.5*(std::sqrt(5.) + 1); // golden ratio
   const double eps = std::pow(10.0, -std::numeric_limits<double>::digits10);

   const std::complex<double> zero(0.0, 0.0);
   CHECK_CLOSE_COMPLEX(Li3(zero), 0.0, eps);

   const std::complex<double> one(1.0, 0.0);
   CHECK_CLOSE_COMPLEX(Li3(one), zeta3, eps);

   const std::complex<double> mone(-1.0, 0.0);
   CHECK_CLOSE_COMPLEX(Li3(mone), -3./4.*zeta3, eps);

   const std::complex<double> half(0.5, 0.0);
   CHECK_CLOSE_COMPLEX(Li3(half), pow3(ln2)/6. - pi2/12.*ln2 + 7./8.*zeta3, eps);

   const std::complex<double> gr(1/sqr(phi), 0.0);
   CHECK_CLOSE_COMPLEX(Li3(gr), 4./5.*zeta3 + 2./3.*pow3(std::log(phi)) - 2./15.*pi2*std::log(phi), eps);
}

// tests signbit for 0.0 and -0.0 arguments
TEST_CASE("test_signed_zero")
{
   // skip test if platform does not supprt signed zero
   if (!has_signed_zero()) {
      return;
   }

   using polylogarithm::Li3;

   const double pz64 = 0.0, nz64 = -0.0;
   const long double pz128 = 0.0L, nz128 = -0.0L;

   // real Li3
   CHECK( std::signbit(Li3(nz64)));
   CHECK(!std::signbit(Li3(pz64)));
   CHECK( std::signbit(poly_Li3(nz64)));
   CHECK(!std::signbit(poly_Li3(pz64)));
#ifdef ENABLE_FORTRAN
   CHECK( std::signbit(poly_Li3_fortran(nz64)));
   CHECK(!std::signbit(poly_Li3_fortran(pz64)));
#endif

   // complex Li3
   CHECK( std::signbit(std::real(Li3(std::complex<double>(nz64, nz64)))));
   CHECK( std::signbit(std::imag(Li3(std::complex<double>(nz64, nz64)))));
   CHECK(!std::signbit(std::real(Li3(std::complex<double>(pz64, nz64)))));
   CHECK( std::signbit(std::imag(Li3(std::complex<double>(pz64, nz64)))));
   CHECK( std::signbit(std::real(Li3(std::complex<double>(nz64, pz64)))));
   CHECK(!std::signbit(std::imag(Li3(std::complex<double>(nz64, pz64)))));
   CHECK(!std::signbit(std::real(Li3(std::complex<double>(pz64, pz64)))));
   CHECK(!std::signbit(std::imag(Li3(std::complex<double>(pz64, pz64)))));

   CHECK( std::signbit(std::real(poly_Li3(std::complex<double>(nz64, nz64)))));
   CHECK( std::signbit(std::imag(poly_Li3(std::complex<double>(nz64, nz64)))));
   CHECK(!std::signbit(std::real(poly_Li3(std::complex<double>(pz64, nz64)))));
   CHECK( std::signbit(std::imag(poly_Li3(std::complex<double>(pz64, nz64)))));
   CHECK( std::signbit(std::real(poly_Li3(std::complex<double>(nz64, pz64)))));
   CHECK(!std::signbit(std::imag(poly_Li3(std::complex<double>(nz64, pz64)))));
   CHECK(!std::signbit(std::real(poly_Li3(std::complex<double>(pz64, pz64)))));
   CHECK(!std::signbit(std::imag(poly_Li3(std::complex<double>(pz64, pz64)))));

#ifdef ENABLE_FORTRAN
   CHECK( std::signbit(std::real(poly_Li3_fortran(std::complex<double>(nz64, nz64)))));
   CHECK( std::signbit(std::imag(poly_Li3_fortran(std::complex<double>(nz64, nz64)))));
   CHECK(!std::signbit(std::real(poly_Li3_fortran(std::complex<double>(pz64, nz64)))));
   CHECK( std::signbit(std::imag(poly_Li3_fortran(std::complex<double>(pz64, nz64)))));
   CHECK( std::signbit(std::real(poly_Li3_fortran(std::complex<double>(nz64, pz64)))));
   CHECK(!std::signbit(std::imag(poly_Li3_fortran(std::complex<double>(nz64, pz64)))));
   CHECK(!std::signbit(std::real(poly_Li3_fortran(std::complex<double>(pz64, pz64)))));
   CHECK(!std::signbit(std::imag(poly_Li3_fortran(std::complex<double>(pz64, pz64)))));
#endif

   CHECK( std::signbit(std::real(Li3(std::complex<long double>(nz128, nz128)))));
   CHECK( std::signbit(std::imag(Li3(std::complex<long double>(nz128, nz128)))));
   CHECK(!std::signbit(std::real(Li3(std::complex<long double>(pz128, nz128)))));
   CHECK( std::signbit(std::imag(Li3(std::complex<long double>(pz128, nz128)))));
   CHECK( std::signbit(std::real(Li3(std::complex<long double>(nz128, pz128)))));
   CHECK(!std::signbit(std::imag(Li3(std::complex<long double>(nz128, pz128)))));
   CHECK(!std::signbit(std::real(Li3(std::complex<long double>(pz128, pz128)))));
   CHECK(!std::signbit(std::imag(Li3(std::complex<long double>(pz128, pz128)))));

   CHECK( std::signbit(std::real(poly_Li3(std::complex<long double>(nz128, nz128)))));
   CHECK( std::signbit(std::imag(poly_Li3(std::complex<long double>(nz128, nz128)))));
   CHECK(!std::signbit(std::real(poly_Li3(std::complex<long double>(pz128, nz128)))));
   CHECK( std::signbit(std::imag(poly_Li3(std::complex<long double>(pz128, nz128)))));
   CHECK( std::signbit(std::real(poly_Li3(std::complex<long double>(nz128, pz128)))));
   CHECK(!std::signbit(std::imag(poly_Li3(std::complex<long double>(nz128, pz128)))));
   CHECK(!std::signbit(std::real(poly_Li3(std::complex<long double>(pz128, pz128)))));
   CHECK(!std::signbit(std::imag(poly_Li3(std::complex<long double>(pz128, pz128)))));
}

template<typename T>
struct Data {
   std::complex<T> z;
   std::complex<T> li_expected;
   T eps;
};

TEST_CASE("test_overflow")
{
   using polylogarithm::Li3;

   const double eps64 = std::pow(10.0, -std::numeric_limits<double>::digits10);
   const long double eps128 = std::pow(10.0L, -std::numeric_limits<long double>::digits10);

   const Data<double> data64[] = {
      {{1e300, 1.0}, {-5.4934049431527088e7, 749538.186928224}, eps64},
      {{1.0, 1e300}, {-5.4936606061973454e7, 374771.031356405}, eps64}
   };

   for (const auto& d : data64) {
      CHECK_CLOSE_COMPLEX(Li3(d.z), d.li_expected, d.eps);
      CHECK_CLOSE_COMPLEX(poly_Li3(d.z), d.li_expected, d.eps);
#ifdef ENABLE_FORTRAN
      CHECK_CLOSE_COMPLEX(poly_Li3_fortran(d.z), d.li_expected, d.eps);
#endif
   }

   const Data<long double> data128[] = {
      {{1e4000L, 1.0L}, {-1.3021939960597718633480560501920418668419e11L, 1.3325123323168432929414697360944689612e8L}, eps128}
   };

   for (const auto& d : data128) {
      CHECK_CLOSE_COMPLEX(Li3(d.z), d.li_expected, d.eps);
      CHECK_CLOSE_COMPLEX(poly_Li3(d.z), d.li_expected, d.eps);
   }
}

TEST_CASE("test_real_fixed_values")
{
   const auto eps64 = std::pow(10.0 , -std::numeric_limits<double>::digits10);

   const std::string filename(std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "Li3.txt");
   const auto fixed_values = polylogarithm::test::read_from_file<long double>(filename);

   for (auto v: fixed_values) {
      const auto z128 = v.first;
      const auto z64 = to<double>(z128);
      const auto x64 = std::real(z64);
      const auto li128_expected = std::real(v.second);
      const auto li64_expected = static_cast<double>(li128_expected);

      if (std::imag(z128) == 0.0L) {
         const auto li64_poly   = polylogarithm::Li3(x64);
         const auto li64_poly_c = li3(x64);
#ifdef ENABLE_FORTRAN
         const auto li64_poly_f = poly_Li3_fortran(x64);
#endif

         INFO("x(64)         = " << x64);
         INFO("Li3(64)  real = " << li64_expected  << " (expected)");
         INFO("Li3(64)  real = " << li64_poly      << " (polylogarithm C++)");
         INFO("Li3(64)  real = " << li64_poly_c    << " (polylogarithm C)");
#ifdef ENABLE_FORTRAN
         INFO("Li3(64)  real = " << li64_poly_f    << " (polylogarithm Fortran)");
#endif

         CHECK_CLOSE(li64_poly  , li64_expected, eps64);
         CHECK_CLOSE(li64_poly_c, li64_expected, eps64);
#ifdef ENABLE_FORTRAN
         CHECK_CLOSE(li64_poly_f, li64_expected, eps64);
#endif
      }
   }
}

TEST_CASE("test_complex_fixed_values")
{
   const auto eps64  = std::pow(10.0 , -std::numeric_limits<double>::digits10);
   const auto eps128 = std::pow(10.0L, -std::numeric_limits<long double>::digits10);

   const std::string filename(std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "Li3.txt");
   const auto fixed_values = polylogarithm::test::read_from_file<long double>(filename);

   for (auto v: fixed_values) {
      const auto z128 = v.first;
      const auto z64 = to<double>(z128);
      const auto li128_expected = v.second;
      const auto li64_expected = to<double>(li128_expected);

      const auto li64_cmpl    = polylogarithm::Li3(z64);
      const auto li64_cmpl_c  = poly_Li3(z64);
#ifdef ENABLE_FORTRAN
      const auto li64_cmpl_f  = poly_Li3_fortran(z64);
#endif
      const auto li128_cmpl   = polylogarithm::Li3(z128);
      const auto li128_cmpl_c = poly_Li3(z128);
      const auto li128_tsil   = tsil_Li3(z128);

      INFO("z(128)        = " << z128);
      INFO("Li3(64)  cmpl = " << li64_expected  << " (expected)");
      INFO("Li3(64)  cmpl = " << li64_cmpl      << " (polylogarithm C++)");
      INFO("Li3(64)  cmpl = " << li64_cmpl_c    << " (polylogarithm C)");
#ifdef ENABLE_FORTRAN
      INFO("Li3(64)  cmpl = " << li64_cmpl_f    << " (polylogarithm Fortran)");
#endif
      INFO("Li3(128) cmpl = " << li128_expected << " (expected)");
      INFO("Li3(128) cmpl = " << li128_tsil     << " (TSIL)");
      INFO("Li3(128) cmpl = " << li128_cmpl     << " (polylogarithm C++)");
      INFO("Li3(128) cmpl = " << li128_cmpl_c   << " (polylogarithm C)");

      CHECK_CLOSE_COMPLEX(li64_cmpl   , li64_expected , 3*eps64);
      CHECK_CLOSE_COMPLEX(li64_cmpl_c , li64_expected , 3*eps64);
#ifdef ENABLE_FORTRAN
      CHECK_CLOSE_COMPLEX(li64_cmpl_f , li64_expected , 3*eps64);
#endif
      CHECK_CLOSE_COMPLEX(li128_cmpl  , li128_expected, 2*eps128);
      CHECK_CLOSE_COMPLEX(li128_cmpl_c, li128_expected, 2*eps128);

      if (std::abs(std::real(z128) - 1.0L) < 1e-5L &&
          std::abs(std::imag(z128)       ) < 1e-5L) {
         // low precision if z is close to (1.0, 0.0)
         // due to log(real(z)) being not veriy precise for real(z) ~ 1
         CHECK_CLOSE_COMPLEX(li128_tsil, li128_expected, 1000*eps128);
      } else if (std::abs(std::real(z128)) < 2000*eps128 &&
                 std::abs(std::imag(z128)) < 2000*eps128) {
         CHECK_CLOSE_COMPLEX(li128_tsil, std::complex<long double>(0.0L, 0.0L), 1e-12L);
      } else {
         CHECK_CLOSE_COMPLEX(li128_tsil, li128_expected, 2*eps128);
      }

      CHECK_SMALL(Relation_1(z64), 1e5*eps64);
      CHECK_SMALL(Relation_2(z64), 1e5*eps64);
      CHECK_SMALL(Relation_3(z64), 1e5*eps64);
   }
}

TEST_CASE("test_complex_random_values")
{
   using polylogarithm::bench::generate_random_complexes;

   const auto eps = std::pow(10.0, -std::numeric_limits<double>::digits10);
   const auto values = generate_random_complexes<double>(10000, -10, 10);

   for (auto v: values) {
      const std::complex<double> li3 = polylogarithm::Li3(v);
      const std::complex<double> li3_c = poly_Li3(v);
#ifdef ENABLE_FORTRAN
      const std::complex<double> li3_f = poly_Li3_fortran(v);
#endif
      const std::complex<double> li3_tsil = to<double>(tsil_Li3(v));

      CHECK_CLOSE_COMPLEX(li3, li3_tsil, 10*eps);
      CHECK_CLOSE_COMPLEX(li3, li3_c   , 10*eps);
#ifdef ENABLE_FORTRAN
      CHECK_CLOSE_COMPLEX(li3, li3_f   , 10*eps);
#endif
   }
}
