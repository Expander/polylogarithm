#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <gsl/gsl_sf_dilog.h>
#include "dilog.hpp"

const std::complex<double> omega(0.5, std::sqrt(3.)/2.);

std::complex<double> values[] = {
   0., 0.5, 1., 1.5,
   -0., -0.5, -1., -1.5,
   -(std::sqrt(5.) - 1)/2.,
   -(std::sqrt(5.) + 1)/2.,
   (std::sqrt(5.) + 1)/2.,
   (std::sqrt(5.) + 3)/2.,
   omega, omega*omega, 1. + omega, 1./(1. + omega)
};

bool is_equal(double a, double b, double eps = std::numeric_limits<double>::epsilon()) {
   return std::abs(a - b) < eps;
}

std::complex<double> clog(std::complex<double> z) {
   std::complex<double> zf(z);
   // convert -0.0 to 0.0
   if (std::real(zf) == 0.0) zf.real() = 0.0;
   if (std::imag(zf) == 0.0) zf.imag() = 0.0;
   return std::log(zf);
}

std::size_t run_real_tests() {
   std::size_t errors = 0;

   for (std::size_t i = 0; i < sizeof(values)/sizeof(values[0]); i++) {
      const double x = std::real(values[i]);
      const double li2 = dilogarithm::dilog(x);
      gsl_sf_result li2_gsl;
      gsl_sf_dilog_e(x, &li2_gsl);

      const bool eq = is_equal(li2, li2_gsl.val, 1e-10);

      if (!eq) {
         std::cout << "Error: dilogarithms are not equal: "
                   << li2 << "(dilog) != "
                   << li2_gsl.val << "(GSL)\n";
         errors++;
      }
   }

   return errors;
}

std::size_t run_complex_tests() {
   std::size_t errors = 0;

   for (std::size_t i = 0; i < sizeof(values)/sizeof(values[0]); i++) {
      const std::complex<double> z = values[i];
      const std::complex<double> li2 = dilogarithm::dilog(z);
      gsl_sf_result li2_gsl_re, li2_gsl_im;
      gsl_sf_complex_dilog_e(std::abs(z), std::arg(z), &li2_gsl_re, &li2_gsl_im);

      const bool eq_re = is_equal(std::real(li2), li2_gsl_re.val, 1e-8);
      const bool eq_im = is_equal(std::imag(li2), li2_gsl_im.val, 1e-8);

      if (!eq_re) {
         std::cout << "Error: real parts of dilogarithms are not equal: "
                   << std::real(li2) << "(dilog) != "
                   << li2_gsl_re.val << "(GSL)\n";
         errors++;
      }

      if (!eq_im) {
         std::cout << "Error: imaginary parts of dilogarithms are not equal: "
                   << std::imag(li2) << "(dilog) != "
                   << li2_gsl_im.val << "(GSL)\n";
         errors++;
      }
   }

   return errors;
}

struct Relation_1 {
   static std::string name() { return "1"; }
   std::complex<double> operator()(std::complex<double> z) const {
      return dilogarithm::dilog(z) + dilogarithm::dilog(-z)
         - dilogarithm::dilog(z*z)/2.;
   }
};

struct Relation_2 {
   static std::string name() { return "2"; }
   std::complex<double> operator()(std::complex<double> z) const {
      return dilogarithm::dilog(1.-z) + dilogarithm::dilog(1.-1./z)
         + clog(z)*clog(z)/2.;
   }
};

struct Relation_3 {
   static std::string name() { return "3"; }
   std::complex<double> operator()(std::complex<double> z) const {
      return dilogarithm::dilog(z) + dilogarithm::dilog(1.-z)
         - (M_PI*M_PI/6. - clog(z) * clog(1.-z));
   }
};

struct Relation_4 {
   static std::string name() { return "4"; }
   std::complex<double> operator()(std::complex<double> z) const {
      return dilogarithm::dilog(-z) - dilogarithm::dilog(1.-z)
         + dilogarithm::dilog(1.-z*z)/2.
         - (-M_PI*M_PI/12. - clog(z) * clog(1.+z));
   }
};

struct Relation_5 {
   static std::string name() { return "5"; }
   std::complex<double> operator()(std::complex<double> z) const {
      return dilogarithm::dilog(z) + dilogarithm::dilog(1./z)
         - (-M_PI*M_PI/6. - clog(-z)*clog(-z)/2.);
   }
};

template <class T>
std::size_t check_relation(std::complex<double> z) {
   const std::complex<double> result(T()(z));

   if (std::real(result) > 1e-9) {
      std::cout << "Error: relation " << T::name()
                << " violated for z = " << z << ": "
                << result << '\n';
      return 1;
   }

   return 0;
}

std::size_t run_relation_tests() {
   std::size_t errors = 0;

   for (std::size_t i = 0; i < sizeof(values)/sizeof(values[0]); i++) {
      const std::complex<double> z = values[i];

      errors += check_relation<Relation_1>(z);
      errors += check_relation<Relation_2>(z);
      errors += check_relation<Relation_3>(z);
      errors += check_relation<Relation_4>(z);
      errors += check_relation<Relation_5>(z);
   }

   return errors;
}

int main() {
   std::size_t errors = 0;

   errors += run_real_tests();
   errors += run_complex_tests();
   errors += run_relation_tests();

   return errors != 0;
}
