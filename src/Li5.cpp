// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the GNU Lesser General Public
// License (GNU LGPL) version 3.
// ====================================================================

#include "Li5.hpp"
#include <cmath>

namespace polylogarithm {

namespace {
   template <typename T>
   T sqr(T x) noexcept { return x*x; }

   template <typename T>
   T pow4(T x) noexcept { return x*x*x*x; }

   // converts -0.0 to 0.0
   std::complex<double> clog(std::complex<double> z) noexcept {
      if (std::real(z) == 0.0) z.real(0.0);
      if (std::imag(z) == 0.0) z.imag(0.0);
      return std::log(z);
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
      0.13953189300411522633744855967078, -0.028633777006172839506172839506173,
      0.0040317412551440329218106995884774, -0.00033985018004115226337448559670782,
      4.5445184621617666490944600374313e-6, 2.3916808048569011882908870275245e-6,
      -1.2762692600122746588544351898175e-7, -3.1628984306505932440256787279547e-8,
      3.2848118445335191621518574238472e-9, 4.7613713995660579048319132897732e-10,
      -8.0846898171909830256460260362332e-11, -7.238764858773720694682921583759e-12,
      1.943976011517396849305564920936e-12, 1.025697840597723597188135594332e-13,
      -4.6180551009884830180582086241066e-14, -1.1535857196470580036842511483419e-15,
      1.0903545401333393987977088380966e-15, 2.3148136317292526394079710319009e-18,
      -2.5669917043265292194334891993397e-17, 4.5708620607314969014495962686014e-19,
      6.0366779613205705882356103311411e-19, -2.167762494406241295879417172184e-20,
      -1.4194096615600165298332266882011e-20, 7.5020009506413862553237752161953e-22,
      3.3387045395078397164371515925447e-22, -2.3060040442620347682521515135259e-23,
      -7.8581732456894818904499064631503e-24, 6.6683453043738808548651370461306e-25,
      1.8509156540925297189464979688365e-25, -1.8591529445174085584103184057636e-26,
      -4.362974648034589044726608174371e-27, 5.0611076099529284482263489587811e-28,
      1.0291918249756878203788897930001e-28, -1.3551391221018316615687785328304e-29,
      -2.4294059612957382655924154095657e-30, 3.5851973966503705211516473805307e-31,
      5.7379658161039720640084653867328e-32, -9.4003593624568734535277452038948e-33
   };

   if (z == 0.)
      return 0.;
   if (z == 1.)
      return zeta5;
   if (z == -1.)
      return -15.*zeta5/16.;

   std::complex<double> u, sum;

   if (std::abs(z) <= 1.) {
      u = -clog(1. - z);
   } else { // az > 1.
      const std::complex<double> lnz  = clog(-z);
      const std::complex<double> lnz2 = sqr(lnz);
      const std::complex<double> lnz4 = pow4(lnz);
      u = -clog(1. - 1./z);
      sum = -1./360.*lnz*(7*PI4 + 10.*PI2*lnz2 + 3.*lnz4);
   }

   std::complex<double> p = 1.;

   for (const double b: bf) {
      p *= u;
      sum += b*p;
   }

   return sum;
}

} // namespace polylogarithm
