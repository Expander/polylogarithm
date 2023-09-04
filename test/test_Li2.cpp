#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "alt.h"
#include "bench.hpp"
#include "c_wrappers.h"
#include "fortran_wrappers.h"
#include "Li2.hpp"
#include "read_data.hpp"
#include <cmath>
#include <complex>
#include <limits>
#include <random>
#include <vector>

#ifndef M_PI
#define M_PI 3.1415926535897932
#endif

#ifdef ENABLE_GSL

#include <gsl/gsl_sf_dilog.h>

double gsl_Li2(double x) {
   gsl_sf_result li2_gsl{};
   gsl_sf_dilog_e(x, &li2_gsl);
   return li2_gsl.val;
}

std::complex<double> gsl_Li2(std::complex<double> z) {
   gsl_sf_result li2_gsl_re{}, li2_gsl_im{};
   gsl_sf_complex_dilog_e(std::abs(z), std::arg(z), &li2_gsl_re, &li2_gsl_im);
   return {li2_gsl_re.val, li2_gsl_im.val};
}

#endif

#define CHECK_CLOSE(a,b,eps) CHECK((a) == doctest::Approx(b).epsilon(eps))
#define CHECK_CLOSE_COMPLEX(a,b,eps) do {                               \
      CHECK_CLOSE(std::real(a), std::real(b), (eps));                   \
      CHECK_CLOSE(std::imag(a), std::imag(b), (eps));                   \
   } while (0)
#define CHECK_SMALL(a,eps) CHECK(std::abs(a) < (eps))

const std::complex<double> omega(0.5, std::sqrt(3.)/2.);
const std::complex<double> zero(0.,0.);

template <class T> T sqr(T x) { return x*x; }

bool is_unity(std::complex<long double> z, long double eps)
{
   return std::abs(std::real(z) - 1.0L) <= eps && std::imag(z) == 0.0L;
}

/// special values to be checked
const std::vector<std::complex<double>> special_values = {
   { 0.0, 0.0 },
   { 0.5, 0.0 },
   { 1.0, 0.0 },
   { 1.5, 0.0 },
   { -0.0, 0.0 },
   { -0.5, 0.0 },
   { -1.0, 0.0 },
   { -1.5, 0.0 },
   { -(std::sqrt(5.) - 1)/2.0, 0.0 },
   { -(std::sqrt(5.) + 1)/2.0, 0.0 },
   { (std::sqrt(5.) + 1)/2.0, 0.0 },
   { (std::sqrt(5.) + 3)/2.0, 0.0 },
   omega,
   omega*omega,
   1. + omega,
   1./(1. + omega)
};

std::complex<double> clog(std::complex<double> z) {
   std::complex<double> zf(z);
   // convert -0.0 to 0.0
   if (std::real(zf) == 0.0) { zf.real(0.0); }
   if (std::imag(zf) == 0.0) { zf.imag(0.0); }
   return std::log(zf);
}

std::complex<float> to_c32(std::complex<long double> z)
{
   return std::complex<float>(std::real(z), std::imag(z));
}

std::complex<double> to_c64(std::complex<long double> z)
{
   return std::complex<double>(std::real(z), std::imag(z));
}

std::complex<double> hdecay_Li2(std::complex<double> z) {
   double li2_re{}, li2_im{};
   hdecay_dilog(std::real(z), std::imag(z), &li2_re, &li2_im);
   return { li2_re, li2_im };
}

std::complex<double> hollik_Li2(std::complex<double> z) {
   double li2_re{}, li2_im{};
   hollik_dilog(std::real(z), std::imag(z), &li2_re, &li2_im);
   return { li2_re, li2_im };
}

std::complex<long double> lt_Li2(std::complex<long double> z) {
   long double re{}, im{};
   feynhiggs_dilog(std::real(z), std::imag(z), &re, &im);
   return { re, im };
}

inline float poly_Li2(float z) {
   return li2f(z);
}

inline double poly_Li2(double z) {
   return li2(z);
}

inline long double poly_Li2(long double z) {
   return li2l(z);
}

std::complex<double> poly_Li2(std::complex<double> z) {
   double re{}, im{};
   cli2_c(std::real(z), std::imag(z), &re, &im);
   return { re, im };
}

#ifdef ENABLE_FORTRAN

double poly_Li2_fortran(double x) {
   double res{};
   li2_fortran(&x, &res);
   return res;
}

std::complex<double> poly_Li2_fortran(std::complex<double> z) {
   const double re = std::real(z);
   const double im = std::imag(z);
   double res_re{}, res_im{};
   cli2_fortran(&re, &im, &res_re, &res_im);
   return { res_re, res_im };
}

#endif

std::complex<long double> poly_Li2(std::complex<long double> z) {
   long double re{}, im{};
   cli2l_c(std::real(z), std::imag(z), &re, &im);
   return { re, im };
}

std::complex<double> sherpa_Li2(std::complex<double> z) {
   double re{}, im{};
   sherpa_dilog(std::real(z), std::imag(z), &re, &im);
   return { re, im };
}

std::complex<double> spheno_Li2(std::complex<double> z) {
   double re{}, im{};
   spheno_dilog(std::real(z), std::imag(z), &re, &im);
   return { re, im };
}

std::complex<long double> tsil_Li2(std::complex<long double> z) {
   long double re{}, im{};
   tsil_dilog_complex(std::real(z), std::imag(z), &re, &im);
   return { re, im };
}

const auto Relation_1 = [](std::complex<double> z) {
   using polylogarithm::Li2;
   return Li2(z) + Li2(-z) - Li2(z*z)/2.;
};

const auto Relation_2 = [](std::complex<double> z) {
   using polylogarithm::Li2;

   if (std::abs(z) < 1e-10 || std::real(z) < 0.) {
      return zero;
   }

   return Li2(1.-z) + Li2(1.-1./z) + clog(z)*clog(z)/2.;
};

const auto Relation_3 = [](std::complex<double> z) {
   using polylogarithm::Li2;

   if (std::abs(z) < 1e-10 || std::abs(std::real(z) - 1.) < 1e-10) {
      return zero;
   }

   return Li2(z) + Li2(1.-z)
      - (M_PI*M_PI/6. - clog(z) * clog(1.-z));
};

const auto Relation_4 = [](std::complex<double> z) {
   using polylogarithm::Li2;

   if (std::abs(z) < 1e-10 || std::abs(std::real(z) + 1.) < 1e-10
       || std::real(z) < 0. || std::imag(z) < 0.) {
      return zero;
   }

   return Li2(-z) - Li2(1.-z) + Li2(1.-z*z)/2.
      - (-M_PI*M_PI/12. - clog(z) * clog(1.+z));
};

const auto Relation_5 = [](std::complex<double> z) {
   using polylogarithm::Li2;

   if (std::abs(z) < 1e-10
       || (std::real(z) > 0. && std::real(z) < 1.)) {
      return zero;
   }

   return Li2(z) + Li2(1./z)
      - (-M_PI*M_PI/6. - clog(-z)*clog(-z)/2.);
};

TEST_CASE("test_special_values")
{
   using polylogarithm::Li2;
   using std::log;
   const double eps = std::pow(10.0 , -std::numeric_limits<double>::digits10);
   const double pi  = M_PI;
   const double pi2 = sqr(pi);
   const double ln2 = std::log(2.);
   const double ln3 = std::log(3.);

   // special values
   CHECK_CLOSE(Li2(-1.), -pi2/12., eps);

   CHECK_CLOSE(Li2(0.), 0., eps);

   CHECK_CLOSE(Li2(0.5), pi2/12. - 0.5*sqr(ln2), eps);

   CHECK_CLOSE(Li2(1.), pi2/6., eps);

   {
      const auto i = std::complex<double>(0.,1.);
      const auto two = std::complex<double>(2.,0.);
      CHECK_CLOSE_COMPLEX(Li2(two), pi2/4. - i*pi*ln2, eps);
   }

   CHECK_CLOSE(Li2(-(sqrt(5.)-1.)/2.),
               -pi2/15. + 0.5*sqr(log((sqrt(5.)-1.)/2.)), eps);

   CHECK_CLOSE(Li2(-(sqrt(5.)+1.)/2.),
               -pi2/10. - sqr(log((sqrt(5.)+1.)/2.)), eps);

   CHECK_CLOSE(Li2((3.-sqrt(5.))/2.),
               pi2/15. - sqr(log((sqrt(5.)-1.)/2.)), eps);

   CHECK_CLOSE(Li2((sqrt(5.)-1.)/2.),
               pi2/10. - sqr(log((sqrt(5.)-1.)/2.)), eps);

   // special value identities
   CHECK_CLOSE(Li2(1./3.) - Li2(1./9.)/6.,
               pi2/18. - sqr(ln3)/6., eps);

   CHECK_CLOSE(Li2(-0.5) + Li2(1./9.)/6.,
               -pi2/18. + ln2*ln3
               - sqr(ln2)/2. - sqr(ln3)/3., eps);

   CHECK_CLOSE(Li2(0.25) + Li2(1./9.)/3.,
               pi2/18. + 2.*ln2*ln3
               - 2.*sqr(ln2) - 2.*sqr(ln3)/3., eps);

   CHECK_CLOSE(Li2(-1./3.) - Li2(1./9.)/3.,
               -pi2/18. + sqr(ln3)/6., eps);

   CHECK_CLOSE(Li2(-1./8.) + Li2(1./9.),
               - 0.5*sqr(log(9./8.)), eps);

   CHECK_CLOSE(36.*Li2(0.5) - 36.*Li2(0.25)
               - 12.*Li2(1./8.) + 6.*Li2(1./64.),
               pi2, 2*eps);

   {
      // special point wher z^2 < epsilon
      const std::complex<double> z(-1.08371e-08,1.32716e-24);
      const std::complex<double> li2_expected(-1.08370999706393160389154078878181e-8, 1.3271599928087172e-24);
      const std::complex<double> li2_poly_cpp(Li2(z));
      const std::complex<double> li2_poly_c(poly_Li2(z));
#ifdef ENABLE_FORTRAN
      const std::complex<double> li2_poly_f(poly_Li2_fortran(z));
#endif

      CHECK_CLOSE_COMPLEX(li2_poly_cpp, li2_expected, eps);
      CHECK_CLOSE_COMPLEX(li2_poly_c  , li2_expected, eps);
#ifdef ENABLE_FORTRAN
      CHECK_CLOSE_COMPLEX(li2_poly_f  , li2_expected, eps);
#endif

      CHECK(std::abs(std::real(li2_poly_cpp - li2_expected)) == 0.0);
      CHECK(std::abs(std::imag(li2_poly_cpp - li2_expected)) == 0.0);
      CHECK(std::abs(std::real(li2_poly_c   - li2_expected)) == 0.0);
      CHECK(std::abs(std::imag(li2_poly_c   - li2_expected)) == 0.0);
#ifdef ENABLE_FORTRAN
      CHECK(std::abs(std::real(li2_poly_f   - li2_expected)) == 0.0);
      CHECK(std::abs(std::imag(li2_poly_f   - li2_expected)) == 0.0);
#endif
   }

   {
      // special point where Re[Li2[z]] == 0
      const std::complex<double> z0(12.5951703698450161286398965, 0.0);
      const std::complex<double> li0(
         -4.41928863051485104279391911740417734355378e-16,
         -7.9586388813396972514450338597465302668858497479134219184303);

      CHECK_CLOSE(Li2(std::real(z0))          , std::real(li0), eps);
      CHECK_CLOSE(algorithm_327(std::real(z0)), std::real(li0), eps);
      CHECK_CLOSE(algorithm_490(std::real(z0)), std::real(li0), eps);
#ifdef ENABLE_GSL
      CHECK_CLOSE(gsl_Li2(std::real(z0))      , std::real(li0), eps);
#endif

      CHECK_CLOSE_COMPLEX(Li2(z0)       , li0, eps);
#ifdef ENABLE_GSL
      CHECK_CLOSE_COMPLEX(gsl_Li2(z0)   , li0, eps);
#endif
      CHECK_CLOSE_COMPLEX(hdecay_Li2(z0), li0, eps);
      CHECK_CLOSE_COMPLEX(hollik_Li2(z0), li0, eps);
      CHECK_CLOSE_COMPLEX(lt_Li2(z0)    , li0, eps);
      CHECK_CLOSE_COMPLEX(tsil_Li2(z0)  , li0, eps);
      CHECK_CLOSE_COMPLEX(sherpa_Li2(z0), li0, 2*eps);
      CHECK_CLOSE_COMPLEX(spheno_Li2(z0), li0, eps);
   }
}

TEST_CASE("test_real_fixed_values")
{
   const auto eps32  = std::pow(10.0f, -std::numeric_limits<float>::digits10);
   const auto eps64  = std::pow(10.0 , -std::numeric_limits<double>::digits10);
   const auto eps128 = std::pow(10.0L, -std::numeric_limits<long double>::digits10);

   const std::string filename(std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "Li2.txt");
   const auto fixed_values = polylogarithm::test::read_from_file<long double>(filename);

   for (auto v: fixed_values) {
      const auto z128 = v.first;
      const auto z32 = to_c32(z128);
      const auto z64 = to_c64(z128);
      const auto x32 = std::real(z32);
      const auto x64 = std::real(z64);
      const auto x128 = std::real(z128);
      const auto li128_expected = std::real(v.second);
      const auto li32_expected = static_cast<float>(li128_expected);
      const auto li64_expected = static_cast<double>(li128_expected);

      if (std::imag(z128) == 0.0L) {
         const auto li64_327      = algorithm_327(x64);
         const auto li64_490      = algorithm_490(x64);
         const auto li64_babar    = babar_dilog(x64);
         const auto li64_cephes   = cephes_dilog(x64);
         const auto li64_cephes_2 = cephes_dilog_2(x64);
         const auto li64_hassani  = hassani_dilog(x64);
         const auto li64_koelbig  = koelbig_dilog(x64);
         const auto li128_koelbig = koelbig_dilogl(x128);
         const auto li64_morris   = morris_dilog(x64);
#ifdef ENABLE_GSL
         const auto li64_gsl      = gsl_Li2(x64);
#endif
         const auto li32_poly     = polylogarithm::Li2(x32);
         const auto li64_poly     = polylogarithm::Li2(x64);
         const auto li128_poly    = polylogarithm::Li2(x128);
         const auto li32_poly_c   = poly_Li2(x32);
         const auto li64_poly_c   = poly_Li2(x64);
         const auto li128_poly_c  = poly_Li2(x128);
#ifdef ENABLE_FORTRAN
         const auto li64_poly_f   = poly_Li2_fortran(x64);
#endif

         INFO("x(32)         = " << x32);
         INFO("Li2(32)  real = " << li32_expected  << " (expected)");
         INFO("Li2(32)  real = " << li32_poly      << " (polylogarithm C++)");
         INFO("Li2(32)  real = " << li32_poly_c    << " (polylogarithm C)");
         INFO("------------------------------------------------------------");
         INFO("x(64)         = " << x64);
         INFO("Li2(64)  real = " << li64_expected  << " (expected)");
         INFO("Li2(64)  real = " << li64_327       << " (algorithm 327)");
         INFO("Li2(64)  real = " << li64_490       << " (algorithm 490)");
         INFO("Li2(64)  real = " << li64_babar     << " (BaBar)");
         INFO("Li2(64)  real = " << li64_cephes    << " (cephes)");
         INFO("Li2(64)  real = " << li64_cephes_2  << " (cephes 2)");
         INFO("Li2(64)  real = " << li64_hassani   << " (hassani)");
         INFO("Li2(64)  real = " << li64_koelbig   << " (koelbig)");
         INFO("Li2(64)  real = " << li64_morris    << " (morris)");
#ifdef ENABLE_GSL
         INFO("Li2(64)  real = " << li64_gsl       << " (GSL)");
#endif
         INFO("Li2(64)  real = " << li64_poly      << " (polylogarithm C++)");
         INFO("Li2(64)  real = " << li64_poly_c    << " (polylogarithm C)");
#ifdef ENABLE_FORTRAN
         INFO("Li2(64)  real = " << li64_poly_f    << " (polylogarithm Fortran)");
#endif
         INFO("------------------------------------------------------------");
         INFO("x(128)        = " << x128);
         INFO("Li2(128) real = " << li128_expected << " (expected)");
         INFO("Li2(128) real = " << li128_poly     << " (polylogarithm C++)");
         INFO("Li2(128) real = " << li128_poly_c   << " (polylogarithm C)");
         INFO("Li2(128) real = " << li128_koelbig  << " (koelbig)");

         CHECK_CLOSE(li32_poly    , std::real(li32_expected) , 2*eps32);
         CHECK_CLOSE(li32_poly_c  , std::real(li32_expected) , 2*eps32);
         CHECK_CLOSE(li64_327     , std::real(li64_expected) , 10*eps64);
         CHECK_CLOSE(li64_490     , std::real(li64_expected) , 2*eps64);
         CHECK_CLOSE(li64_babar   , std::real(li64_expected) , 100*eps64);
         CHECK_CLOSE(li64_cephes  , std::real(li64_expected) , 2*eps64);
         CHECK_CLOSE(li64_cephes_2, std::real(li64_expected) , 2*eps64);
         CHECK_CLOSE(li64_hassani , std::real(li64_expected) , 100*eps64);
         CHECK_CLOSE(li64_koelbig , std::real(li64_expected) , 2*eps64);
         CHECK_CLOSE(li64_morris  , std::real(li64_expected) , 2*eps64);
#ifdef ENABLE_GSL
         CHECK_CLOSE(li64_gsl     , std::real(li64_expected) , 2*eps64);
#endif
         CHECK_CLOSE(li64_poly    , std::real(li64_expected) , 2*eps64);
         CHECK_CLOSE(li64_poly_c  , std::real(li64_expected) , 2*eps64);
#ifdef ENABLE_FORTRAN
         CHECK_CLOSE(li64_poly_f  , std::real(li64_expected) , 2*eps64);
#endif
         CHECK_CLOSE(li128_poly   , std::real(li128_expected), 2*eps128);
         CHECK_CLOSE(li128_poly_c , std::real(li128_expected), 2*eps128);
         CHECK_CLOSE(li128_koelbig, std::real(li128_expected), 2*eps128);
      }
   }
}

TEST_CASE("test_complex_fixed_values")
{
   const auto eps32  = std::pow(10.0f, -std::numeric_limits<float>::digits10);
   const auto eps64  = std::pow(10.0 , -std::numeric_limits<double>::digits10);
   const auto eps128 = std::pow(10.0L, -std::numeric_limits<long double>::digits10);

   const std::string filename(std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "Li2.txt");
   const auto fixed_values = polylogarithm::test::read_from_file<long double>(filename);

   for (auto v: fixed_values) {
      const auto z128 = v.first;
      const auto z32 = to_c32(z128);
      const auto z64 = to_c64(z128);
      const auto li128_expected = v.second;
      const auto li32_expected = to_c32(li128_expected);
      const auto li64_expected = to_c64(li128_expected);

      const auto li32_poly   = polylogarithm::Li2(z32);
      const auto li64_poly   = polylogarithm::Li2(z64);
      const auto li128_poly  = polylogarithm::Li2(z128);
      const auto li64_poly_c = poly_Li2(z64);
      const auto li128_poly_c= poly_Li2(z128);
#ifdef ENABLE_FORTRAN
      const auto li64_poly_f = poly_Li2_fortran(z64);
#endif
#ifdef ENABLE_GSL
      const auto li64_gsl    = gsl_Li2(z64);
#endif
      const auto li64_hdecay = hdecay_Li2(z64);
      const auto li64_hollik = hollik_Li2(z64);
      const auto li128_lt    = lt_Li2(z128);
      const auto li128_tsil  = tsil_Li2(z128);
      const auto li64_sherpa = sherpa_Li2(z64);
      const auto li64_spheno = spheno_Li2(z64);

      INFO("z(32)         = " << z32);
      INFO("Li2(32)  cmpl = " << li32_expected  << " (expected)");
      INFO("Li2(32)  cmpl = " << li32_poly      << " (polylogarithm C++)");
      INFO("------------------------------------------------------------");
      INFO("z(64)         = " << z64);
      INFO("Li2(64)  cmpl = " << li64_expected  << " (expected)");
#ifdef ENABLE_GSL
      INFO("Li2(64)  cmpl = " << li64_gsl       << " (GSL)");
#endif
      INFO("Li2(64)  cmpl = " << li64_hdecay    << " (HDECAY)");
      INFO("Li2(64)  cmpl = " << li64_hollik    << " (Hollik)");
      INFO("Li2(64)  cmpl = " << li64_poly      << " (polylogarithm C++)");
      INFO("Li2(64)  cmpl = " << li64_poly_c    << " (polylogarithm C)");
#ifdef ENABLE_FORTRAN
      INFO("Li2(64)  cmpl = " << li64_poly_f    << " (polylogarithm Fortran)");
#endif
      INFO("Li2(64)  cmpl = " << li64_sherpa    << " (Sherpa)");
      INFO("Li2(64)  cmpl = " << li64_spheno    << " (SPheno)");
      INFO("------------------------------------------------------------");
      INFO("z(128)        = " << z128);
      INFO("Li2(128) cmpl = " << li128_expected << " (expected)");
      INFO("Li2(128) cmpl = " << li128_poly     << " (polylogarithm C++)");
      INFO("Li2(128) cmpl = " << li128_poly_c   << " (polylogarithm C)");
      INFO("Li2(128) cmpl = " << li128_lt       << " (LoopTools)");
      INFO("Li2(128) cmpl = " << li128_tsil     << " (TSIL)");

      if (!(std::real(z32) == 1 && std::imag(z32) < 1e-23)) {
         CHECK_CLOSE_COMPLEX(li32_poly   , li32_expected , 2*eps32);
      }
      CHECK_CLOSE_COMPLEX(li64_poly   , li64_expected , 2*eps64);
      CHECK_CLOSE_COMPLEX(li64_poly_c , li64_expected , 2*eps64);
#ifdef ENABLE_FORTRAN
      CHECK_CLOSE_COMPLEX(li64_poly_f , li64_expected , 2*eps64);
#endif
#ifdef ENABLE_GSL
      CHECK_CLOSE_COMPLEX(li64_gsl    , li64_expected , 10*eps64);
#endif
      CHECK_CLOSE_COMPLEX(li64_hdecay , li64_expected , 2*eps64);
      CHECK_CLOSE_COMPLEX(li64_hollik , li64_expected , 2*eps64);
      CHECK_CLOSE_COMPLEX(li64_sherpa , li64_expected , 2*eps64);
      CHECK_CLOSE_COMPLEX(li64_spheno , li64_expected , 2*eps64);
      CHECK_CLOSE_COMPLEX(li128_poly  , li128_expected, 2*eps128);
      CHECK_CLOSE_COMPLEX(li128_poly_c, li128_expected, 2*eps128);
      CHECK_CLOSE_COMPLEX(li128_lt    , li128_expected, 2*eps128);
      CHECK_CLOSE_COMPLEX(li128_tsil  , li128_expected, 4*eps128);
   }
}

TEST_CASE("test_real_random_values")
{
   using polylogarithm::bench::generate_random_scalars;

   const auto eps64  = std::pow(10.0 , -std::numeric_limits<double>::digits10);
   const auto values = generate_random_scalars<double>(10000, -10, 10);

   for (auto v: values) {
      const double li2   = polylogarithm::Li2(v);
      const double li2_c = poly_Li2(v);
#ifdef ENABLE_FORTRAN
      const double li2_f = poly_Li2(v);
#endif
#ifdef ENABLE_GSL
      const double li2_gsl = gsl_Li2(v);
#endif
      const double li2_327 = algorithm_327(v);
      const double li2_490 = algorithm_490(v);
      const double li2_babar = babar_dilog(v);
      const double li2_cephes = cephes_dilog(v);
      const double li2_cephes_2 = cephes_dilog_2(v);
      const double li2_hassani = hassani_dilog(v);
      const double li2_koelbig = koelbig_dilog(v);
      const double li2_morris = morris_dilog(v);

      INFO("x = " << v);
      INFO("Li2(64) real = " << li2          << " (polylogarithm C++)");
      INFO("Li2(64) real = " << li2_c        << " (polylogarithm C)");
#ifdef ENABLE_FORTRAN
      INFO("Li2(64) real = " << li2_f        << " (polylogarithm Fortran)");
#endif
#ifdef ENABLE_GSL
      INFO("Li2(64) real = " << li2_gsl      << " (GSL)");
#endif
      INFO("Li2(64) real = " << li2_327      << " (algorithm 327)");
      INFO("Li2(64) real = " << li2_490      << " (algorithm 490)");
      INFO("Li2(64) real = " << li2_babar    << " (BaBar)");
      INFO("Li2(64) real = " << li2_cephes   << " (cephes)");
      INFO("Li2(64) real = " << li2_cephes_2 << " (cephes 2)");
      INFO("Li2(64) real = " << li2_hassani  << " (Hassani)");
      INFO("Li2(64) real = " << li2_koelbig  << " (Koelbig)");
      INFO("Li2(64) real = " << li2_morris   << " (Morris)");

      CHECK_CLOSE(li2, li2_c       , eps64);
#ifdef ENABLE_FORTRAN
      CHECK_CLOSE(li2, li2_f       , eps64);
#endif
#ifdef ENABLE_GSL
      CHECK_CLOSE(li2, li2_gsl     , 2*eps64);
#endif
      CHECK_CLOSE(li2, li2_327     , 10*eps64);
      CHECK_CLOSE(li2, li2_490     , eps64);
      CHECK_CLOSE(li2, li2_babar   , 100*eps64);
      CHECK_CLOSE(li2, li2_cephes  , 2*eps64);
      CHECK_CLOSE(li2, li2_cephes_2, 2*eps64);
      CHECK_CLOSE(li2, li2_hassani , 100*eps64);
      CHECK_CLOSE(li2, li2_koelbig , 2*eps64);
      CHECK_CLOSE(li2, li2_morris  , eps64);
   }
}

TEST_CASE("test_complex_random_values")
{
   using polylogarithm::bench::generate_random_complexes;

   const auto eps = std::pow(10.0, -std::numeric_limits<double>::digits10);
   const auto values = generate_random_complexes<double>(10000, -10, 10);

   for (auto v: values) {
      const std::complex<double> li2   = polylogarithm::Li2(v);
      const std::complex<double> li2_c = poly_Li2(v);
#ifdef ENABLE_FORTRAN
      const std::complex<double> li2_f = poly_Li2_fortran(v);
#endif
#ifdef ENABLE_GSL
      const std::complex<double> li2_gsl = gsl_Li2(v);
#endif
      const std::complex<double> li2_hdecay = hdecay_Li2(v);
      const std::complex<double> li2_hollik = hollik_Li2(v);
      const std::complex<double> li2_sherpa = sherpa_Li2(v);
      const std::complex<double> li2_spheno = spheno_Li2(v);
      const std::complex<double> li2_tsil = to_c64(tsil_Li2(v));
      const std::complex<double> li2_lt   = to_c64(lt_Li2(v));

      INFO("z = " << v);
      INFO("Li2(64) cmpl = " << li2        << " (polylogarithm C++)");
      INFO("Li2(64) cmpl = " << li2_c      << " (polylogarithm C)");
#ifdef ENABLE_FORTRAN
      INFO("Li2(64) cmpl = " << li2_f      << " (polylogarithm Fortran)");
#endif
#ifdef ENABLE_GSL
      INFO("Li2(64) cmpl = " << li2_gsl    << " (GSL)");
#endif
      INFO("Li2(64) cmpl = " << li2_hdecay << " (HDECAY)");
      INFO("Li2(64) cmpl = " << li2_hollik << " (Hollik)");
      INFO("Li2(64) cmpl = " << li2_lt     << " (LoopTools)");
      INFO("Li2(64) cmpl = " << li2_sherpa << " (Sherpa)");
      INFO("Li2(64) cmpl = " << li2_spheno << " (SPheno)");
      INFO("Li2(64) cmpl = " << li2_tsil   << " (TSIL)");

#ifdef ENABLE_GSL
      CHECK_CLOSE_COMPLEX(li2, li2_gsl   , 10*eps);
#endif
      CHECK_CLOSE_COMPLEX(li2, li2_c     , 10*eps);
#ifdef ENABLE_FORTRAN
      CHECK_CLOSE_COMPLEX(li2, li2_f     , 10*eps);
#endif
      CHECK_CLOSE_COMPLEX(li2, li2_tsil  , 10*eps);
      CHECK_CLOSE_COMPLEX(li2, li2_hdecay, 10*eps);
      CHECK_CLOSE_COMPLEX(li2, li2_hollik, 10*eps);
      CHECK_CLOSE_COMPLEX(li2, li2_lt    , 2*eps);
      CHECK_CLOSE_COMPLEX(li2, li2_sherpa, 10*eps);
      CHECK_CLOSE_COMPLEX(li2, li2_spheno, 10*eps);
   }
}

TEST_CASE("test_relations")
{
   const double eps = 2*std::pow(10.0, -std::numeric_limits<double>::digits10);

   for (const auto v: special_values) {
      CHECK_SMALL(Relation_1(v), eps);
      CHECK_SMALL(Relation_2(v), eps);
      CHECK_SMALL(Relation_3(v), eps);
      CHECK_SMALL(Relation_4(v), eps);
      CHECK_SMALL(Relation_5(v), eps);
   }
}
