#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "dilog.hpp"
#include <cmath>
#include <gsl/gsl_sf_dilog.h>
#include <iostream>
#include <limits>
#include <random>
#include <string>

#define CHECK_CLOSE(a,b,eps) CHECK((a) == doctest::Approx(b).epsilon(eps))
#define CHECK_SMALL(a,eps) CHECK(std::abs(a) < (eps))

const std::complex<double> omega(0.5, std::sqrt(3.)/2.);
const std::complex<double> zero(0.,0.);

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

double gsl_dilog(double x) {
   gsl_sf_result li2_gsl{};
   gsl_sf_dilog_e(x, &li2_gsl);
   return li2_gsl.val;
}

std::complex<double> gsl_dilog(std::complex<double> z) {
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
   return dilogarithm::dilog(z) + dilogarithm::dilog(-z)
      - dilogarithm::dilog(z*z)/2.;
};

const auto Relation_2 = [](std::complex<double> z) {
   if (std::abs(z) < 1e-10 || std::real(z) < 0.)
      return zero;

   return dilogarithm::dilog(1.-z) + dilogarithm::dilog(1.-1./z)
      + clog(z)*clog(z)/2.;
};

const auto Relation_3 = [](std::complex<double> z) {
   if (std::abs(z) < 1e-10 || std::abs(std::real(z) - 1.) < 1e-10)
      return zero;

   return dilogarithm::dilog(z) + dilogarithm::dilog(1.-z)
      - (M_PI*M_PI/6. - clog(z) * clog(1.-z));
};

const auto Relation_4 = [](std::complex<double> z) {
   if (std::abs(z) < 1e-10 || std::abs(std::real(z) + 1.) < 1e-10
       || std::real(z) < 0. || std::imag(z) < 0.)
      return zero;

   return dilogarithm::dilog(-z) - dilogarithm::dilog(1.-z)
      + dilogarithm::dilog(1.-z*z)/2.
      - (-M_PI*M_PI/12. - clog(z) * clog(1.+z));
};

const auto Relation_5 = [](std::complex<double> z) {
   if (std::abs(z) < 1e-10
       || (std::real(z) > 0. && std::real(z) < 1.))
      return zero;

   return dilogarithm::dilog(z) + dilogarithm::dilog(1./z)
      - (-M_PI*M_PI/6. - clog(-z)*clog(-z)/2.);
};

template <class T>
std::size_t check_relation(T rel, std::complex<double> z)
{
   const std::complex<double> result = rel(z);
   CHECK_SMALL(std::real(result), 1e-9);
   CHECK_SMALL(std::imag(result), 1e-9);
}

TEST_CASE("test_real_fixed_values")
{
   for (auto v: values) {
      const double x = std::real(v);
      const double li2 = dilogarithm::dilog(x);
      const double li2_gsl = gsl_dilog(x);

      CHECK_CLOSE(li2, li2_gsl, 1e-10);
   }
}

TEST_CASE("test_real_random_values")
{
   const auto values = generate_random_doubles(10000, -10, 10);

   for (auto v: values) {
      const double x = std::real(v);
      const double li2 = dilogarithm::dilog(x);
      const double li2_gsl = gsl_dilog(x);

      CHECK_CLOSE(li2, li2_gsl, 1e-10);
   }
}

TEST_CASE("test_complex_fixed_values")
{
   for (auto v: values) {
      const std::complex<double> li2 = dilogarithm::dilog(v);
      const std::complex<double> li2_gsl = gsl_dilog(v);

      CHECK_CLOSE(std::real(li2), std::real(li2_gsl), 1e-8);
      CHECK_CLOSE(std::imag(li2), std::imag(li2_gsl), 1e-8);
   }
}

TEST_CASE("test_complex_random_values")
{
   const auto values = generate_random_complexes(10000, -10, 10);

   for (auto v: values) {
      const std::complex<double> li2 = dilogarithm::dilog(v);
      const std::complex<double> li2_gsl = gsl_dilog(v);

      CHECK_CLOSE(std::real(li2), std::real(li2_gsl), 1e-8);
      CHECK_CLOSE(std::imag(li2), std::imag(li2_gsl), 1e-8);
   }
}

TEST_CASE("test_relations")
{
   for (auto v: values) {
      const std::complex<double> z = v;

      check_relation(Relation_1, z);
      check_relation(Relation_2, z);
      check_relation(Relation_3, z);
      check_relation(Relation_4, z);
      check_relation(Relation_5, z);
   }
}
