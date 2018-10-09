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
   const double PI = 3.1415926535897932384626433832795;
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
   const double PI    = 3.1415926535897932384626433832795;
   const double PI2   = PI*PI;
   const double PI4   = PI2*PI2;
   const double zeta5 = 1.0369277551433699263313654864570;
   static const int N = 40;

   const double bf[N] = {
      1., -15./32.,
      0.13953189300411522633744855967078   , -0.028633777006172839506172839506173  ,
      0.0040317412551440329218106995884774 , -0.00033985018004115226337448559670782,
      4.5445184621617666490944600374313e-6 ,  2.3916808048569011882908870275245e-6 ,
     -1.2762692600122746588544351898175e-7 , -3.1628984306505932440256787279547e-8 ,
      3.2848118445335191621518574238472e-9 ,  4.7613713995660579048319132897732e-10,
     -8.0846898171909830256460260362332e-11, -7.238764858773720694682921583759e-12 ,
      1.943976011517396849305564920936e-12 ,  1.025697840597723597188135594332e-13 ,
     -4.6180551009884830180582086241066e-14, -1.1535857196470580036842511483419e-15,
      1.0903545401333393987977088380966e-15,  2.3148136317292526394079710319009e-18,
     -2.5669917043265292194334891993397e-17,  4.5708620607314969014495962686014e-19,
      6.0366779613205705882356103311411e-19, -2.167762494406241295879417172184e-20 ,
     -1.4194096615600165298332266882011e-20,  7.5020009506413862553237752161953e-22,
      3.3387045395078397164371515925447e-22, -2.3060040442620347682521515135259e-23,
     -7.8581732456894818904499064631503e-24,  6.6683453043738808548651370461306e-25,
      1.8509156540925297189464979688365e-25, -1.8591529445174085584103184057636e-26,
     -4.362974648034589044726608174371e-27 ,  5.0611076099529284482263489587811e-28,
      1.0291918249756878203788897930001e-28, -1.3551391221018316615687785328304e-29,
     -2.4294059612957382655924154095657e-30,  3.5851973966503705211516473805307e-31,
      5.7379658161039720640084653867328e-32, -9.4003593624568734535277452038948e-33
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
      const double zeta3 = 1.2020569031595942853997381615114;
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
