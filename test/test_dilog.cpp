#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "dilog.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <gsl/gsl_sf_dilog.h>

#define CHECK_CLOSE(a,b,eps) CHECK((a) == doctest::Approx(b).epsilon(eps))
#define CHECK_SMALL(a,eps) CHECK((a) < eps)

const std::complex<double> omega(0.5, std::sqrt(3.)/2.);

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

const auto Relation_1 = [](std::complex<double> z) {
   return dilogarithm::dilog(z) + dilogarithm::dilog(-z)
      - dilogarithm::dilog(z*z)/2.;
};

const auto Relation_2 = [](std::complex<double> z) {
   if (std::abs(z) < 1e-10)
      return std::complex<double>(0.,0.);

   return dilogarithm::dilog(1.-z) + dilogarithm::dilog(1.-1./z)
      + clog(z)*clog(z)/2.;
};

const auto Relation_3 = [](std::complex<double> z) {
   if (std::abs(z) < 1e-10 || std::abs(std::real(z) - 1.) < 1e-10)
      return std::complex<double>(0.,0.);

   return dilogarithm::dilog(z) + dilogarithm::dilog(1.-z)
      - (M_PI*M_PI/6. - clog(z) * clog(1.-z));
};

const auto Relation_4 = [](std::complex<double> z) {
   if (std::abs(z) < 1e-10 || std::abs(std::real(z) + 1.) < 1e-10)
      return std::complex<double>(0.,0.);

   return dilogarithm::dilog(-z) - dilogarithm::dilog(1.-z)
      + dilogarithm::dilog(1.-z*z)/2.
      - (-M_PI*M_PI/12. - clog(z) * clog(1.+z));
};

const auto Relation_5 = [](std::complex<double> z) {
   if (std::abs(z) < 1e-10)
      return std::complex<double>(0.,0.);

   return dilogarithm::dilog(z) + dilogarithm::dilog(1./z)
      - (-M_PI*M_PI/6. - clog(-z)*clog(-z)/2.);
};

template <class T>
std::size_t check_relation(T r, std::complex<double> z)
{
   const std::complex<double> result = r(z);
   CHECK_SMALL(std::real(result), 1e-9);
}

TEST_CASE("test_real")
{
   for (auto v: values) {
      const double x = std::real(v);
      const double li2 = dilogarithm::dilog(x);
      gsl_sf_result li2_gsl;
      gsl_sf_dilog_e(x, &li2_gsl);

      CHECK_CLOSE(li2, li2_gsl.val, 1e-10);
   }
}

TEST_CASE("test_complex")
{
   for (auto v: values) {
      const std::complex<double> z = v;
      const std::complex<double> li2 = dilogarithm::dilog(z);
      gsl_sf_result li2_gsl_re, li2_gsl_im;
      gsl_sf_complex_dilog_e(std::abs(z), std::arg(z), &li2_gsl_re, &li2_gsl_im);

      CHECK_CLOSE(std::real(li2), li2_gsl_re.val, 1e-8);
      CHECK_CLOSE(std::imag(li2), li2_gsl_im.val, 1e-8);
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
