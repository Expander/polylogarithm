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
   template <typename T> T pow4(T x) noexcept { return x*x*x*x; }

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
 * @brief Clausen function \f$\mathrm{Cl}_5(\theta) = \mathrm{Re}(\mathrm{Li}_5(e^{i\theta}))\f$
 * @param x real angle
 * @return \f$\mathrm{Cl}_5(\theta)\f$
 */
double Cl5(double x)
{
   using std::exp;
   const double PI = 3.141592653589793;
   const std::complex<double> i(0.,1.);

   while (x >= 2*PI)
      x -= 2*PI;

   while (x < 0.)
      x += 2*PI;

   return std::real(Li5(exp(i*x)));
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
   static const int N = 40;

   const double bf[N] = {
      1., -15./32.,
      1.395318930041152e-01, -2.863377700617283e-02,
      4.031741255144032e-03, -3.398501800411522e-04,
      4.544518462161766e-06,  2.391680804856901e-06,
     -1.276269260012274e-07, -3.162898430650593e-08,
      3.284811844533519e-09,  4.761371399566057e-10,
     -8.084689817190983e-11, -7.238764858773720e-12,
      1.943976011517396e-12,  1.025697840597723e-13,
     -4.618055100988483e-14, -1.153585719647058e-15,
      1.090354540133339e-15,  2.314813631729252e-18,
     -2.566991704326529e-17,  4.570862060731496e-19,
      6.036677961320570e-19, -2.167762494406241e-20,
     -1.419409661560016e-20,  7.502000950641386e-22,
      3.338704539507839e-22, -2.306004044262034e-23,
     -7.858173245689481e-24,  6.668345304373880e-25,
      1.850915654092529e-25, -1.859152944517408e-26,
     -4.362974648034589e-27,  5.061107609952928e-28,
      1.029191824975687e-28, -1.355139122101831e-29,
     -2.429405961295738e-30,  3.585197396650370e-31,
      5.737965816103972e-32, -9.400359362456873e-33
   };

   if (is_close(z, 0.))
      return 0.;
   if (is_close(z, 1.))
      return zeta5;
   if (is_close(z, 1., 0.02)) {
      const std::complex<double> I(0.,1.);
      const std::complex<double> IPI(0.,PI);
      const std::complex<double> zm1 = z - 1.;
      const std::complex<double> lzm1 = clog(zm1);
      const double zeta3 = 1.202056903159594;
      const double ceil = std::arg(zm1) > 0. ? 1. : 0.;
      const std::complex<double> cs[8] = {
         PI4/90.,
         (-PI4/180. + zeta3/2.),
         (15.*PI2 + 2.*PI4 - 270.*zeta3)/540.,
         (-lzm1/24. + (125. + 60.*(-1 + 2*ceil)*IPI - 4.*PI2*(15 + PI2) + 660.*zeta3)/1440.),
         (lzm1/12. + (-565. + 300.*(1 - 2*ceil)*IPI + 175.*PI2 + 8.*PI4 - 1500.*zeta3)/3600.),
         (-17./144.*lzm1 + (1779. + 1020.*(-1 + 2*ceil)*IPI - 450.*PI2 - 16.*PI4 + 3288.*zeta3)/8640.),
         (-0.23923611111111112 + 7./48.*(1 - 2*ceil)*IPI + 7./48.*lzm1 + 29./540.*PI2
          + PI4/630. - 7./20.*zeta3),
         (-967./5760.*lzm1 + (1266861. + 812280.*(-1 + 2*ceil)*IPI
                              - 560.*PI2*(469 + 12*PI2) + 1568160.*zeta3)/4.8384e6)
      };

      return zeta5 +
         zm1 * (cs[0] +
         zm1 * (cs[1] +
         zm1 * (cs[2] +
         zm1 * (cs[3] +
         zm1 * (cs[4] +
         zm1 * (cs[5] +
         zm1 * (cs[6] +
         zm1 * (cs[7]))))))));
   }
   if (is_close(z, -1.))
      return -15.*zeta5/16.;

   std::complex<double> u, rest;

   if (std::abs(z) <= 1.) {
      u = -clog(1. - z);
   } else { // az > 1.
      const std::complex<double> lnz  = clog(-z);
      const std::complex<double> lnz2 = pow2(lnz);
      const std::complex<double> lnz4 = pow4(lnz);
      u = -clog(1. - 1./z);
      rest = -1./360.*lnz*(7*PI4 + 10.*PI2*lnz2 + 3.*lnz4);
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
      u * (bf[18] +
      u * (bf[19] +
      u * (bf[20] +
      u * (bf[21] +
      u * (bf[22] +
      u * (bf[23] +
      u * (bf[24] +
      u * (bf[25] +
      u * (bf[26] +
      u * (bf[27] +
      u * (bf[28] +
      u * (bf[29] +
      u * (bf[30] +
      u * (bf[31] +
      u * (bf[32] +
      u * (bf[33] +
      u * (bf[34] +
      u * (bf[35] +
      u * (bf[36] +
      u * (bf[37] +
      u * (bf[38] +
      u * (bf[39]))))))))))))))))))))))))))))))))))))))));
}

} // namespace polylogarithm
