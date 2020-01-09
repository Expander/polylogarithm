#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "Li2.hpp"
#include "bench.hpp"
#include "tsil_cpp.h"
#include "read_data.hpp"
#include <cmath>
#include <complex>
#include <gsl/gsl_sf_dilog.h>
#include <random>
#include <vector>

#define CHECK_CLOSE(a,b,eps) CHECK((a) == doctest::Approx(b).epsilon(eps))
#define CHECK_CLOSE_COMPLEX(a,b,eps) do {                               \
      CHECK_CLOSE(std::real(a), std::real(b), (eps));                   \
      CHECK_CLOSE(std::imag(a), std::imag(b), (eps));                   \
   } while (0)
#define CHECK_SMALL(a,eps) CHECK(std::abs(a) < (eps))

const std::complex<double> omega(0.5, std::sqrt(3.)/2.);
const std::complex<double> zero(0.,0.);

template <class T> T sqr(T x) { return x*x; }

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
   if (std::real(zf) == 0.0) zf.real(0.0);
   if (std::imag(zf) == 0.0) zf.imag(0.0);
   return std::log(zf);
}

std::complex<double> to_c64(std::complex<long double> z)
{
   return std::complex<double>(std::real(z), std::imag(z));
}

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

std::complex<long double> tsil_Li2(std::complex<long double> z) {
   return TSIL_Dilog_(z);
}

const auto Relation_1 = [](std::complex<double> z) {
   using polylogarithm::Li2;
   return Li2(z) + Li2(-z) - Li2(z*z)/2.;
};

const auto Relation_2 = [](std::complex<double> z) {
   using polylogarithm::Li2;

   if (std::abs(z) < 1e-10 || std::real(z) < 0.)
      return zero;

   return Li2(1.-z) + Li2(1.-1./z) + clog(z)*clog(z)/2.;
};

const auto Relation_3 = [](std::complex<double> z) {
   using polylogarithm::Li2;

   if (std::abs(z) < 1e-10 || std::abs(std::real(z) - 1.) < 1e-10)
      return zero;

   return Li2(z) + Li2(1.-z)
      - (M_PI*M_PI/6. - clog(z) * clog(1.-z));
};

const auto Relation_4 = [](std::complex<double> z) {
   using polylogarithm::Li2;

   if (std::abs(z) < 1e-10 || std::abs(std::real(z) + 1.) < 1e-10
       || std::real(z) < 0. || std::imag(z) < 0.)
      return zero;

   return Li2(-z) - Li2(1.-z) + Li2(1.-z*z)/2.
      - (-M_PI*M_PI/12. - clog(z) * clog(1.+z));
};

const auto Relation_5 = [](std::complex<double> z) {
   using polylogarithm::Li2;

   if (std::abs(z) < 1e-10
       || (std::real(z) > 0. && std::real(z) < 1.))
      return zero;

   return Li2(z) + Li2(1./z)
      - (-M_PI*M_PI/6. - clog(-z)*clog(-z)/2.);
};

TEST_CASE("test_special_values")
{
   using polylogarithm::Li2;
   using std::log;
   const double eps = 1e-15;
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
               pi2, eps);
}

TEST_CASE("test_fixed_values")
{
   const auto eps64  = 1e-15;
   const auto eps128 = 1e-15L; // @todo increase precision
   const std::string filename(std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "Li2.txt");
   const auto fixed_values = polylogarithm::test::read_from_file<long double>(filename);

   for (auto v: fixed_values) {
      const auto z128 = v.first;
      const auto z64 = to_c64(z128);
      const auto x64 = std::real(z64);
      const auto x128 = std::real(z128);
      const auto li128_expected = v.second;
      const auto li64_expected = to_c64(li128_expected);

      const auto li64_real  = polylogarithm::Li2(x64);
      const auto li128_real = polylogarithm::Li2(x128);
      const auto li64_cmpl  = polylogarithm::Li2(z64);
      const auto li128_cmpl = polylogarithm::Li2(z128);
      const auto li64_gsl   = gsl_Li2(x64);
      const auto li128_tsil = tsil_Li2(z128);

      INFO("z(128)        = " << z128);
      INFO("Li2(64)  cmpl = " << li64_expected << " (expected)");
      INFO("Li2(128) cmpl = " << li128_expected << " (expected)");
      INFO("Li2(128) cmpl = " << li128_tsil << " (TSIL)");
      INFO("Li2(64)  cmpl = " << li64_cmpl << " (polylogarithm)");
      INFO("Li2(128) cmpl = " << li128_cmpl << " (polylogarithm)");

      CHECK_CLOSE_COMPLEX(li64_cmpl , li64_expected , eps64);
      CHECK_CLOSE_COMPLEX(li128_tsil, li128_expected, eps128);
      CHECK_CLOSE_COMPLEX(li128_cmpl, li128_expected, eps128);

      if (std::imag(z128) == 0.0L) {
         INFO("Li2(64)  real = " << li64_gsl << " (GSL)");
         INFO("Li2(64)  real = " << li64_real << " (polylogarithm)");
         INFO("Li2(128) real = " << li128_real << " (polylogarithm)");

         CHECK_CLOSE(li64_real , std::real(li64_expected) , eps64);
         CHECK_CLOSE(li64_gsl  , std::real(li64_expected) , eps64);
         CHECK_CLOSE(li128_real, std::real(li128_expected), eps128);
      }
   }
}

TEST_CASE("test_real_random_values")
{
   using namespace polylogarithm::bench;

   const auto values = generate_random_doubles(10000, -10, 10);

   for (auto v: values) {
      const double li2 = polylogarithm::Li2(v);
      const double li2_gsl = gsl_Li2(v);

      CHECK_CLOSE(li2, li2_gsl, 1e-15);
   }
}

TEST_CASE("test_complex_random_values")
{
   using namespace polylogarithm::bench;

   const auto values = generate_random_complexes(10000, -10, 10);

   for (auto v: values) {
      const std::complex<double> li2 = polylogarithm::Li2(v);
      const std::complex<double> li2_gsl = gsl_Li2(v);
      const std::complex<double> li2_tsil = to_c64(tsil_Li2(v));

      CHECK_CLOSE_COMPLEX(li2, li2_gsl, 1e-13);
      CHECK_CLOSE_COMPLEX(li2, li2_tsil, 1e-13);
   }
}

TEST_CASE("test_relations")
{
   for (const auto v: special_values) {
      CHECK_SMALL(Relation_1(v), 1e-15);
      CHECK_SMALL(Relation_2(v), 1e-15);
      CHECK_SMALL(Relation_3(v), 1e-15);
      CHECK_SMALL(Relation_4(v), 1e-15);
      CHECK_SMALL(Relation_5(v), 1e-15);
   }
}
