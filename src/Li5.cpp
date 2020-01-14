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
   const double eps_d = 10.0*std::numeric_limits<double>::epsilon();

   template <typename T> T pow2(T x) noexcept { return x*x; }

   // converts -0.0 to 0.0
   std::complex<double> clog(std::complex<double> z) noexcept {
      if (std::real(z) == 0.0) { z.real(0.0); }
      if (std::imag(z) == 0.0) { z.imag(0.0); }
      return std::log(z);
   }

   bool is_close(const std::complex<double>& a, double b, double eps)
   {
      return std::abs(std::real(a) - b) < eps && std::abs(std::imag(a)) < eps;
   }

   template <typename T>
   std::complex<T> cadd(T a, const std::complex<T>& b) noexcept
   {
      return std::complex<T>(a + std::real(b), std::imag(b));
   }

   template <typename T>
   std::complex<T> cadd(const std::complex<T>& a, const std::complex<T>& b) noexcept
   {
      return std::complex<T>(std::real(a) + std::real(b),
                             std::imag(a) + std::imag(b));
   }

   template <typename T>
   std::complex<T> cmul(const std::complex<T>& a, T b) noexcept
   {
      return std::complex<T>(std::real(a) * b, std::imag(a) * b);
   }

   template <typename T>
   std::complex<T> cmul(const std::complex<T>& a, const std::complex<T>& b) noexcept
   {
      return std::complex<T>(
         std::real(a) * std::real(b) - std::imag(a) * std::imag(b),
         std::real(a) * std::imag(b) + std::imag(a) * std::real(b));
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
   const std::complex<double> i(0.0, 1.0);

   while (x >= 2*PI) {
      x -= 2*PI;
   }

   while (x < 0.0) {
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
      1.0, -15.0/32.0,
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

   if (is_close(z, 0.0, eps_d)) {
      return { 0.0, 0.0 };
   }
   if (is_close(z, 1.0, eps_d)) {
      return { zeta5, 0.0 };
   }
   if (is_close(z, -1.0, eps_d)) {
      return { -15.0*zeta5/16.0, 0.0 };
   }

   const auto az  = std::abs(z);
   const auto pz  = std::arg(z);
   const auto lnz = std::log(az);

   if (pow2(lnz) + pow2(pz) < 1.0) { // |log(z)| < 1
      const auto u  = std::complex<double>(lnz, pz); // clog(z)
      const auto u2 = u*u;
      const auto c0 = zeta5;
      const auto c1 = 1.082323233711138; // zeta(4)
      const auto c2 = 0.6010284515797971; // zeta(3)/2
      const auto c3 = 0.2741556778080377;
      const auto c4 = (25.0/12.0 - clog(-u))/24.0;
      const auto c5 = -1.0/240.0;

      const double cs[6] = {
         -1.157407407407407e-04, 2.066798941798942e-07,
         -1.093544413650234e-09, 8.698648744945041e-12,
         -8.689958786158882e-14, 1.008125408021881e-15
      };

      return
         cadd(c0,
         cadd(cmul(u, c1),
         cmul(u2, cadd(c2,
         cadd(cmul(u, c3),
         cmul(u2, cadd(c4,
         cadd(cmul(u, c5),
         cmul(u2, cadd(cs[0],
         cmul(u2, cadd(cs[1],
         cmul(u2, cadd(cs[2],
         cmul(u2, cadd(cs[3],
         cmul(u2, cadd(cs[4],
         cmul(u2, cs[5])))))))))))))))))));
   }

   std::complex<double> u(0.0, 0.0), rest(0.0, 0.0);

   if (az <= 1.0) {
      u = -clog(1.0 - z);
   } else { // az > 1
      auto arg = PI + pz;
      if (arg > PI) { arg -= 2*PI; }
      const auto lmz = std::complex<double>(lnz, arg); // clog(-z)
      const auto lmz2 = pow2(lmz);
      u = -clog(1.0 - 1.0/z);
      rest = -1.0/360.0*lmz*(7*PI4 + lmz2*(10.0*PI2 + 3.0*lmz2));
   }

   return
      cadd(rest,
      cmul(u, cadd(bf[0],
      cmul(u, cadd(bf[1],
      cmul(u, cadd(bf[2],
      cmul(u, cadd(bf[3],
      cmul(u, cadd(bf[4],
      cmul(u, cadd(bf[5],
      cmul(u, cadd(bf[6],
      cmul(u, cadd(bf[7],
      cmul(u, cadd(bf[8],
      cmul(u, cadd(bf[9],
      cmul(u, cadd(bf[10],
      cmul(u, cadd(bf[11],
      cmul(u, cadd(bf[12],
      cmul(u, cadd(bf[13],
      cmul(u, cadd(bf[14],
      cmul(u, cadd(bf[15],
      cmul(u, cadd(bf[16],
      cmul(u, cadd(bf[17],
      cmul(u, bf[18]))))))))))))))))))))))))))))))))))))));
}

} // namespace polylogarithm
