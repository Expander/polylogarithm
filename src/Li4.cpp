// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the GNU Lesser General Public
// License (GNU LGPL) version 3.
// ====================================================================

#include "Li4.hpp"
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
 * @brief Clausen function \f$\mathrm{Cl}_4(\theta) = \mathrm{Im}(\mathrm{Li}_4(e^{i\theta}))\f$
 * @param x real angle
 * @return \f$\mathrm{Cl}_4(\theta)\f$
 */
double Cl4(double x)
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

   return std::imag(Li4(exp(i*x)));
}

/**
 * @brief Complex polylogarithm \f$\mathrm{Li}_4(z)\f$
 * @param z complex argument
 * @return \f$\mathrm{Li}_4(z)\f$
 */
std::complex<double> Li4(const std::complex<double>& z)
{
   const double PI    = 3.141592653589793;
   const double PI2   = PI*PI;
   const double PI4   = PI2*PI2;
   const double zeta4 = 1.082323233711138;
   const double bf[18] = {
      1., -7./16.,
      1.165123456790123e-01, -1.982060185185185e-02,
      1.927932098765432e-03, -3.105709876543209e-05,
     -1.562400911485783e-05,  8.485123546773206e-07,
      2.290961660318971e-07, -2.183261421852691e-08,
     -3.882824879172015e-09,  5.446292103220332e-10,
      6.960805210682725e-11, -1.337573768644521e-11,
     -1.278485268526657e-12,  3.260562858024892e-13,
      2.364757116861825e-14, -7.923135122031161e-15,
   };

   if (is_close(z, 0.))
      return 0.;
   if (is_close(z, 1.))
      return zeta4;
   if (is_close(z, -1.))
      return -7.*PI4/720.;

   const auto az  = std::abs(z);
   const auto pz  = std::arg(z);
   const auto lnz = std::log(az);

   if (pow2(lnz) + pow2(pz) < 1.) { // |log(z)| < 1
      const auto u  = clog(z);
      const auto u2 = u*u;
      const auto c1 = 1.202056903159594; // zeta(3)
      const auto c2 = 0.8224670334241132;
      const auto c3 = (11.0/6.0 - clog(-u))/6.0;
      const auto c4 = -1.0/48.0;

      const double cs[7] = {
         -6.944444444444444e-04, 1.653439153439153e-06,
         -1.093544413650234e-08, 1.043837849393405e-10,
         -1.216594230062244e-12, 1.613000652835010e-14,
         -2.342881045287934e-16
      };

      return zeta4 + u2 * (c2 + u2 * c4) +
         u * (
            c1 +
            u2 * (c3 +
            u2 * (cs[0] +
            u2 * (cs[1] +
            u2 * (cs[2] +
            u2 * (cs[3] +
            u2 * (cs[4] +
            u2 * (cs[5] +
            u2 * (cs[6]))))))))
         );
   }

   std::complex<double> u, r;
   double sgn = 1;

   if (az <= 1.) {
      u = -clog(1. - z);
   } else { // az > 1.
      const auto lnz  = clog(-z);
      const auto lnz2 = pow2(lnz);
      u = -clog(1. - 1./z);
      r = 1./360.*(-7*PI4 + lnz2*(-30.*PI2 - 15.*lnz2));
      sgn = -1;
   }

   const auto sum =
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
