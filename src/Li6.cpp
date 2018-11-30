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
   const double epsilon = std::pow(10., -std::floor(std::numeric_limits<double>::digits10));

   template <typename T> T pow2(T x) noexcept { return x*x; }

   // converts -0.0 to 0.0
   std::complex<double> clog(std::complex<double> z) noexcept {
      if (std::real(z) == 0.0) z.real(0.0);
      if (std::imag(z) == 0.0) z.imag(0.0);
      return std::log(z);
   }

   bool is_close(const std::complex<double>& a, const std::complex<double>& b,
                 double eps = epsilon)
   {
      return std::abs(std::real(a) - std::real(b)) < eps &&
             std::abs(std::imag(a) - std::imag(b)) < eps;
   }
} // anonymous namespace

/**
 * @brief Clausen function \f$\mathrm{Cl}_6(\theta) = \mathrm{Im}(\mathrm{Li}_6(e^{i\theta}))\f$
 * @param x real angle
 * @return \f$\mathrm{Cl}_6(\theta)\f$
 */
double Cl6(double x)
{
   using std::exp;
   const double PI = 3.141592653589793;
   const std::complex<double> i(0.,1.);

   while (x >= 2*PI)
      x -= 2*PI;

   while (x < 0.)
      x += 2*PI;

   if (std::abs(x) < std::numeric_limits<double>::epsilon() ||
       std::abs(x - PI) < std::numeric_limits<double>::epsilon() ||
       std::abs(x - 2*PI) < std::numeric_limits<double>::epsilon())
      return 0.;

   return std::imag(Li6(exp(i*x)));
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
      1., -31./64.,
      1.524134087791495e-01, -3.436555587705761e-02,
      5.717479723936899e-03, -6.818045374657064e-04,
      4.996036194873449e-05, -4.916605119603904e-07,
     -3.063297516130216e-07,  1.441459927084909e-08,
      3.727243823092410e-09, -3.730086734548760e-10,
     -5.124652681608583e-11,  9.054193095663668e-12,
      6.738188261551251e-13, -2.121583115030313e-13,
     -6.840881171901169e-15,  4.869117846200558e-15
   };

   if (is_close(z, 0.))
      return 0.;
   if (is_close(z, 1.))
      return zeta6;
   if (is_close(z, -1.))
      return -31.*zeta6/32.;

   const auto az  = std::abs(z);
   const auto pz  = std::arg(z);
   const auto lnz = std::log(az);

   if (pow2(lnz) + pow2(pz) < 1.) { // |log(z)| < 1
      const auto u  = clog(z);
      const auto u2 = u*u;
      const auto c0 = zeta6;
      const auto c1 = 1.036927755143370; // zeta(5)
      const auto c2 = 0.5411616168555691;
      const auto c3 = 0.2003428171932657;
      const auto c4 = 0.06853891945200943;
      const auto c5 = (137./60. - clog(-u))/120.;
      const auto c6 = -1./1440.;

      const double cs[5] = {
         -1.653439153439153e-05, 2.296443268665491e-08,
         -9.941312851365761e-11, 6.691268265342339e-13,
         -5.793305857439255e-15
      };

      return c0 + u * c1 +
         u2 * (c2 + u * c3 +
         u2 * (c4 + u * c5 +
         u2 * (c6 +
         u * (cs[0] +
         u2 * (cs[1] +
         u2 * (cs[2] +
         u2 * (cs[3] +
         u2 * (cs[4]))))))));
   }

   std::complex<double> u, r;
   double sgn = 1;

   if (std::abs(z) <= 1.) {
      u = -clog(1. - z);
   } else { // az > 1.
      const std::complex<double> lnz  = clog(-z);
      const std::complex<double> lnz2 = pow2(lnz);
      const std::complex<double> lnz4 = pow2(lnz2);
      const std::complex<double> lnz6 = lnz2*lnz4;
      u = -clog(1. - 1./z);
      r = -31.*PI6/15120. - 7./720.*PI4*lnz2 - 1./144.*PI2*lnz4 - 1./720.*lnz6;
      sgn = -1;
   }

   const std::complex<double> sum =
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
      u * (bf[17]))))))))))))))))));

   return sgn*sum + r;
}

} // namespace polylogarithm
