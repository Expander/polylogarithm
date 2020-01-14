// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the GNU Lesser General Public
// License (GNU LGPL) version 3.
// ====================================================================

#include "Li6.hpp"
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
 * @brief Clausen function \f$\mathrm{Cl}_6(\theta) = \mathrm{Im}(\mathrm{Li}_6(e^{i\theta}))\f$
 * @param x real angle
 * @return \f$\mathrm{Cl}_6(\theta)\f$
 */
double Cl6(double x)
{
   const double PI = 3.141592653589793;
   const std::complex<double> i(0.0, 1.0);

   while (x >= 2*PI) {
      x -= 2*PI;
   }

   while (x < 0.0) {
      x += 2*PI;
   }

   if (std::abs(x) < std::numeric_limits<double>::epsilon() ||
       std::abs(x - PI) < std::numeric_limits<double>::epsilon() ||
       std::abs(x - 2*PI) < std::numeric_limits<double>::epsilon()) {
      return 0.0;
   }

   return std::imag(Li6(std::exp(i*x)));
}

/**
 * @brief Complex polylogarithm \f$\mathrm{Li}_6(z)\f$
 * @param z complex argument
 * @return \f$\mathrm{Li}_6(z)\f$
 */
std::complex<double> Li6(const std::complex<double>& z)
{
   const double PI    = 3.141592653589793;
   const double PI2   = PI*PI;
   const double PI4   = PI2*PI2;
   const double PI6   = PI2*PI4;
   const double zeta6 = 1.017343061984449;
   const double bf[18] = {
      1.0, -31.0/64.0,
      1.524134087791495e-01, -3.436555587705761e-02,
      5.717479723936899e-03, -6.818045374657064e-04,
      4.996036194873449e-05, -4.916605119603904e-07,
     -3.063297516130216e-07,  1.441459927084909e-08,
      3.727243823092410e-09, -3.730086734548760e-10,
     -5.124652681608583e-11,  9.054193095663668e-12,
      6.738188261551251e-13, -2.121583115030313e-13,
     -6.840881171901169e-15,  4.869117846200558e-15
   };

   if (is_close(z, 0.0, eps_d)) {
      return { 0.0, 0.0 };
   }
   if (is_close(z, 1.0, eps_d)) {
      return { zeta6, 0.0 };
   }
   if (is_close(z, -1.0, eps_d)) {
      return { -31.0*zeta6/32.0, 0.0 };
   }

   const auto az  = std::abs(z);
   const auto pz  = std::arg(z);
   const auto lnz = std::log(az);

   if (pow2(lnz) + pow2(pz) < 1.0) { // |log(z)| < 1
      const auto u  = std::complex<double>(lnz, pz); // clog(z)
      const auto u2 = u*u;
      const auto c0 = zeta6;
      const auto c1 = 1.036927755143370; // zeta(5)
      const auto c2 = 0.5411616168555691;
      const auto c3 = 0.2003428171932657;
      const auto c4 = 0.06853891945200943;
      const auto c5 = (137.0/60.0 - clog(-u))/120.0;
      const auto c6 = -1.0/1440.0;

      const double cs[5] = {
         -1.653439153439153e-05, 2.296443268665491e-08,
         -9.941312851365761e-11, 6.691268265342339e-13,
         -5.793305857439255e-15
      };

      return
         cadd(c0,
         cadd(cmul(u, c1),
         cmul(u2, cadd(c2,
         cadd(cmul(u, c3),
         cmul(u2, cadd(c4,
         cadd(cmul(u, c5),
         cmul(u2, cadd(c6,
         cmul(u,  cadd(cs[0],
         cmul(u2, cadd(cs[1],
         cmul(u2, cadd(cs[2],
         cmul(u2, cadd(cs[3],
         cmul(u2, cs[4])))))))))))))))))));
   }

   std::complex<double> u(0.0, 0.0), r(0.0, 0.0);
   double sgn = 1;

   if (az <= 1.0) {
      u = -clog(1.0 - z);
   } else { // az > 1
      auto arg = PI + pz;
      if (arg > PI) { arg -= 2*PI; }
      const auto lmz = std::complex<double>(lnz, arg); // clog(-z)
      const auto lmz2 = pow2(lmz);
      u = -clog(1.0 - 1.0/z);
      r = -31.0*PI6/15120.0 + lmz2*(-7.0/720.0*PI4 + lmz2*(-1.0/144.0*PI2 - 1.0/720.0*lmz2));
      sgn = -1;
   }

   const auto sum =
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
      cmul(u, bf[17])))))))))))))))))))))))))))))))))));

   return sgn*sum + r;
}

} // namespace polylogarithm
