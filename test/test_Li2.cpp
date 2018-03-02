#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "Li2.hpp"
#include <cmath>
#include <complex>
#include <gsl/gsl_sf_dilog.h>
#include <random>
#include <vector>

#define CHECK_CLOSE(a,b,eps) CHECK((a) == doctest::Approx(b).epsilon(eps))
#define CHECK_CLOSE_COMPLEX(a,b,eps) do {                               \
      CHECK(std::real(a) == doctest::Approx(std::real(b)).epsilon(eps)); \
      CHECK(std::imag(a) == doctest::Approx(std::imag(b)).epsilon(eps)); \
   } while (0);
#define CHECK_SMALL(a,eps) CHECK(std::abs(a) < (eps))

const std::complex<double> omega(0.5, std::sqrt(3.)/2.);
const std::complex<double> zero(0.,0.);

template <class T> T sqr(T x) { return x*x; }

/// special values to be checked
const std::vector<std::complex<double>> values = {
   0., 0.5, 1., 1.5,
   -0., -0.5, -1., -1.5,
   -(std::sqrt(5.) - 1)/2.,
   -(std::sqrt(5.) + 1)/2.,
   (std::sqrt(5.) + 1)/2.,
   (std::sqrt(5.) + 3)/2.,
   omega, omega*omega, 1. + omega, 1./(1. + omega)
};

std::complex<double> clog(std::complex<double> z) {
   std::complex<double> zf(z);
   // convert -0.0 to 0.0
   if (std::real(zf) == 0.0) zf.real(0.0);
   if (std::imag(zf) == 0.0) zf.imag(0.0);
   return std::log(zf);
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

std::vector<double> generate_random_doubles(int n, double start, double stop)
{
   std::minstd_rand gen;
   std::uniform_real_distribution<double> dist(start, stop);

   std::vector<double> v(n);
   std::generate(begin(v), end(v),
                 [&dist,&gen](){ return dist(gen); });

   return v;
}

std::vector<std::complex<double>> generate_random_complexes(
   int n, double start, double stop)
{
   const auto reals = generate_random_doubles(n, start, stop);
   const auto imags = generate_random_doubles(n, start, stop);

   std::vector<std::complex<double>> v(n);

   for (int i = 0; i < n; i++)
      v[i] = std::complex<double>(reals[i], imags[i]);

   return v;
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

TEST_CASE("test_real_fixed_values")
{
   for (auto v: values) {
      const double x = std::real(v);
      const double li2 = polylogarithm::Li2(x);
      const double li2_gsl = gsl_Li2(x);

      CHECK_CLOSE(li2, li2_gsl, 1e-10);
   }
}

TEST_CASE("test_special_values")
{
   using polylogarithm::Li2;
   using std::log;
   const double eps = 1e-10;
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
      CHECK_CLOSE_COMPLEX(Li2({2.,0.}), pi2/4. - i*pi*ln2, eps);
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

TEST_CASE("test_real_random_values")
{
   const auto values = generate_random_doubles(10000, -10, 10);

   for (auto v: values) {
      const double li2 = polylogarithm::Li2(v);
      const double li2_gsl = gsl_Li2(v);

      CHECK_CLOSE(li2, li2_gsl, 1e-10);
   }
}

TEST_CASE("test_complex_fixed_values")
{
   for (auto v: values) {
      const std::complex<double> li2 = polylogarithm::Li2(v);
      const std::complex<double> li2_gsl = gsl_Li2(v);

      CHECK_CLOSE(std::real(li2), std::real(li2_gsl), 1e-8);
      CHECK_CLOSE(std::imag(li2), std::imag(li2_gsl), 1e-8);
   }
}

TEST_CASE("test_complex_random_values")
{
   const auto values = generate_random_complexes(10000, -10, 10);

   for (auto v: values) {
      const std::complex<double> li2 = polylogarithm::Li2(v);
      const std::complex<double> li2_gsl = gsl_Li2(v);

      CHECK_CLOSE(std::real(li2), std::real(li2_gsl), 1e-8);
      CHECK_CLOSE(std::imag(li2), std::imag(li2_gsl), 1e-8);
   }
}

TEST_CASE("test_relations")
{
   for (const auto v: values) {
      CHECK_SMALL(Relation_1(v), 1e-9);
      CHECK_SMALL(Relation_2(v), 1e-9);
      CHECK_SMALL(Relation_3(v), 1e-9);
      CHECK_SMALL(Relation_4(v), 1e-9);
      CHECK_SMALL(Relation_5(v), 1e-9);
   }
}
