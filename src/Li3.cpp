// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the GNU Lesser General Public
// License (GNU LGPL) version 3.
// ====================================================================

#include "Li3.hpp"
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
 * @brief Clausen function \f$\mathrm{Cl}_3(\theta) = \mathrm{Re}(\mathrm{Li}_3(e^{i\theta}))\f$
 * @param x real angle
 * @return \f$\mathrm{Cl}_3(\theta)\f$
 */
double Cl3(double x)
{
   using std::exp;
   const double PI = 3.1415926535897932384626433832795;
   const std::complex<double> i(0.,1.);

   while (x >= 2*PI)
      x -= 2*PI;

   while (x < 0.)
      x += 2*PI;

   return std::real(Li3(exp(i*x)));
}

/**
 * @brief Complex trilogarithm \f$\mathrm{Li}_3(z)\f$
 * @param z complex argument
 * @return \f$\mathrm{Li}_3(z)\f$
 */
std::complex<double> Li3(const std::complex<double>& z)
{
   const double PI    = 3.141592653589793;
   const double PI2   = PI*PI;
   const double zeta2 = 1.644934066848226;
   const double zeta3 = 1.202056903159594;
   const double inv_e = 0.3678794411714423; // 1/e

   const double bf[40] = {
      1., -3./8., 17./216., -5./576.,
      0.00012962962962962962962962962962963,  0.000081018518518518518518518518518519,
     -3.4193571608537594932152755282007e-6 , -1.3286564625850340136054421768707e-6,
      8.6608717561098513479465860418241e-8 ,  2.5260875955320399764844209288654e-8,
     -2.1446944683640647609338850757365e-9 , -5.14011062201297891533581769272e-10,
      5.2495821146008294363940888085581e-11,  1.0887754406636318375372971570425e-11,
     -1.2779396094493695305581831754072e-12, -2.3698241773087452099797778810124e-13,
      3.1043578879654622942847532704656e-14,  5.2617586299125060841318392511225e-15,
     -7.5384795499492653659925014322677e-16, -1.1862322577752285253082500951246e-16,
      1.8316979965491383382089273121282e-17,  2.7068171031837350151490734712617e-18,
     -4.4554338978296388264326309921763e-19, -6.2375484922556946503653222473984e-20,
      1.0851521534874534913136560996864e-20,  1.4491174866036081930734904966528e-21,
     -2.6466339754458990334740891186144e-22, -3.3897653488510104721925816586081e-23,
      6.4640477336033108890325309821953e-24,  7.975834489602412424209222725905e-25,
     -1.5809178790287483355921117629383e-25, -1.8861499729622868193110225398853e-26,
      3.8715536638418473303997127188831e-27,  4.4801175002345607304865389832051e-28,
     -9.493033871911836126417536760277e-29 , -1.0682813809077381224018214303381e-29,
      2.3304478936103051860078519901928e-30,  2.556077572651975408056356982867e-31,
     -5.7274216061372596844727445803306e-32, -6.1347132137964235825854929689777e-33
   };

   if (is_close(z, 0.))
      return 0.;
   if (is_close(z, 1.))
      return zeta3;
   if (is_close(z, 1., inv_e)) {
      const auto u  = clog(z);
      const auto u2 = u*u;
      const auto u3 = u*u2;
      const auto c0 = zeta3 + zeta2*u - u3/12.;
      const auto c1 = 0.25 * (3.0 - 2.0*clog(-u));

      const std::complex<double> cs[7] = {
         -3.472222222222222e-03, 1.157407407407407e-05,
         -9.841899722852104e-08, 1.148221634332745e-09,
         -1.581572499080917e-11, 2.419500979252515e-13,
         -3.982897776989488e-15
      };

      return c0 +
         u2 * (c1 +
         u2 * (cs[0] +
         u2 * (cs[1] +
         u2 * (cs[2] +
         u2 * (cs[3] +
         u2 * (cs[4] +
         u2 * (cs[5] +
         u2 * (cs[6]))))))));
   }
   if (is_close(z, -1.))
      return -0.75*zeta3;
   if (is_close(z, 0.5)) {
      const double ln2  = 0.6931471805599453; // ln(2)
      const double ln23 = 0.3330246519889295; // ln(2)^3
      return (-2*PI2*ln2 + 4*ln23 + 21*zeta3)/24.;
   }

   std::complex<double> u, rest;

   if (std::abs(z) <= 1.) {
      u = -clog(1. - z);
   } else { // az > 1.
      const auto lmz = clog(-z);
      u = -clog(1. - 1./z);
      rest = -lmz*(pow2(lmz)/6. + PI2/6.);
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
