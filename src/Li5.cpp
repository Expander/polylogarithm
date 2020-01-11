// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the GNU Lesser General Public
// License (GNU LGPL) version 3.
// ====================================================================

#include "Li5.hpp"
#include <cmath>
#include <limits>

namespace polylogarithm {

namespace {
   const double epsilon = std::pow(10., -std::floor(std::numeric_limits<double>::digits10));

   template <typename T> T pow2(T x) noexcept { return x*x; }

   // converts -0.0 to 0.0
   std::complex<double> clog(std::complex<double> z) noexcept {
      if (std::real(z) == 0.0) { z.real(0.0); }
      if (std::imag(z) == 0.0) { z.imag(0.0); }
      return std::log(z);
   }

   bool is_close(const std::complex<double>& a, double b, double eps = epsilon)
   {
      return std::abs(std::real(a) - b) < eps && std::abs(std::imag(a)) < eps;
   }
   } // anonymous namespace

/**
 * @brief Clausen function \f$\mathrm{Cl}_5(\theta) = \mathrm{Re}(\mathrm{Li}_5(e^{i\theta}))\f$
 * @param x real angle
 * @return \f$\mathrm{Cl}_5(\theta)\f$
 */
double Cl5(double x)
{
   const double PI = 3.141592653589793;
   const std::complex<double> i(0.,1.);

   while (x >= 2*PI) {
      x -= 2*PI;
   }

   while (x < 0.) {
      x += 2*PI;
   }

   return std::real(Li5(std::exp(i*x)));
}

/**
 * @brief Complex polylogarithm \f$\mathrm{Li}_5(z)\f$
 * @param z complex argument
 * @return \f$\mathrm{Li}_5(z)\f$
 */
std::complex<double> Li5(const std::complex<double>& z)
{
   const double PI    = 3.141592653589793;
   const double PI2   = PI*PI;
   const double PI4   = PI2*PI2;
   const double zeta5 = 1.036927755143369;
   const double bf[19] = {
      1., -15./32.,
      1.395318930041152e-01, -2.863377700617283e-02,
      4.031741255144032e-03, -3.398501800411522e-04,
      4.544518462161766e-06,  2.391680804856901e-06,
     -1.276269260012274e-07, -3.162898430650593e-08,
      3.284811844533519e-09,  4.761371399566057e-10,
     -8.084689817190983e-11, -7.238764858773720e-12,
      1.943976011517396e-12,  1.025697840597723e-13,
     -4.618055100988483e-14, -1.153585719647058e-15,
      1.090354540133339e-15
   };

   if (is_close(z, 0.)) {
      return { 0., 0. };
   }
   if (is_close(z, 1.)) {
      return { zeta5, 0. };
   }
   if (is_close(z, -1.)) {
      return { -15.*zeta5/16., 0. };
   }

   const auto az  = std::abs(z);
   const auto pz  = std::arg(z);
   const auto lnz = std::log(az);

   if (pow2(lnz) + pow2(pz) < 1.) { // |log(z)| < 1
      const auto u  = std::complex<double>(lnz, pz); // clog(z)
      const auto u2 = u*u;
      const auto c0 = zeta5;
      const auto c1 = 1.082323233711138; // zeta(4)
      const auto c2 = 0.6010284515797971; // zeta(3)/2
      const auto c3 = 0.2741556778080377;
      const auto c4 = (25./12. - clog(-u))/24.;
      const auto c5 = -1./240.;

      const double cs[6] = {
         -1.157407407407407e-04, 2.066798941798942e-07,
         -1.093544413650234e-09, 8.698648744945041e-12,
         -8.689958786158882e-14, 1.008125408021881e-15
      };

      return c0 + u * c1 +
         u2 * (c2 + u * c3 +
         u2 * (c4 + u * c5 +
         u2 * (cs[0] +
         u2 * (cs[1] +
         u2 * (cs[2] +
         u2 * (cs[3] +
         u2 * (cs[4] +
         u2 * (cs[5]))))))));
   }

   std::complex<double> u(0.,0.), rest(0.,0.);

   if (az <= 1.) {
      u = -clog(1. - z);
   } else { // az > 1.
      auto arg = PI + pz;
      if (arg > PI) { arg -= 2*PI; }
      const auto lmz = std::complex<double>(lnz, arg); // clog(-z)
      const auto lmz2 = pow2(lmz);
      u = -clog(1. - 1./z);
      rest = -1./360.*lmz*(7*PI4 + lmz2*(10.*PI2 + 3.*lmz2));
   }

   return rest +
      u * (bf[0] +
      u * (bf[1] +
      u * (bf[2] +
      u * (bf[3] +
      u * (bf[4] +
      u * (bf[5] +
      u * (bf[6] +
      u * (bf[7] +
      u * (bf[8] +
      u * (bf[9] +
      u * (bf[10] +
      u * (bf[11] +
      u * (bf[12] +
      u * (bf[13] +
      u * (bf[14] +
      u * (bf[15] +
      u * (bf[16] +
      u * (bf[17] +
      u * (bf[18])))))))))))))))))));
}

} // namespace polylogarithm
