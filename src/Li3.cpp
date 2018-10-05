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
   template <typename T> T pow3(T x) noexcept { return x*x*x; }
   template <typename T> T pow4(T x) noexcept { return x*x*x*x; }
   template <typename T> T pow5(T x) noexcept { return x*x*x*x*x; }
   template <typename T> T pow6(T x) noexcept { return x*x*x*x*x*x; }
   template <typename T> T pow7(T x) noexcept { return x*x*x*x*x*x*x; }
   template <typename T> T pow8(T x) noexcept { return x*x*x*x*x*x*x*x; }

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
   const double PI    = 3.1415926535897932384626433832795;
   const double PI2   = PI*PI;
   const double zeta3 = 1.2020569031595942853997381615114;
   static const int N = 40;

   const double bf[N] = {
      1., -3./8., 17./216., -5./576.,
      0.00012962962962962962962962962962963, 0.000081018518518518518518518518518519,
      -3.4193571608537594932152755282007e-6, -1.3286564625850340136054421768707e-6,
      8.6608717561098513479465860418241e-8, 2.5260875955320399764844209288654e-8,
      -2.1446944683640647609338850757365e-9, -5.14011062201297891533581769272e-10,
      5.2495821146008294363940888085581e-11, 1.0887754406636318375372971570425e-11,
      -1.2779396094493695305581831754072e-12, -2.3698241773087452099797778810124e-13,
      3.1043578879654622942847532704656e-14, 5.2617586299125060841318392511225e-15,
      -7.5384795499492653659925014322677e-16, -1.1862322577752285253082500951246e-16,
      1.8316979965491383382089273121282e-17, 2.7068171031837350151490734712617e-18,
      -4.4554338978296388264326309921763e-19, -6.2375484922556946503653222473984e-20,
      1.0851521534874534913136560996864e-20, 1.4491174866036081930734904966528e-21,
      -2.6466339754458990334740891186144e-22, -3.3897653488510104721925816586081e-23,
      6.4640477336033108890325309821953e-24, 7.975834489602412424209222725905e-25,
      -1.5809178790287483355921117629383e-25, -1.8861499729622868193110225398853e-26,
      3.8715536638418473303997127188831e-27, 4.4801175002345607304865389832051e-28,
      -9.493033871911836126417536760277e-29, -1.0682813809077381224018214303381e-29,
      2.3304478936103051860078519901928e-30, 2.556077572651975408056356982867e-31,
      -5.7274216061372596844727445803306e-32, -6.1347132137964235825854929689777e-33
   };

   if (is_close(z, 0.))
      return 0.;
   if (is_close(z, 1.))
      return zeta3;
   if (is_close(z, 1., 0.02)) {
      const std::complex<double> I(0.,1.);
      const std::complex<double> IPI(0.,PI);
      const std::complex<double> zm1 = z - 1.;
      const std::complex<double> lzm1 = clog(zm1);
      const double ceil = std::arg(zm1) > 0. ? 1. : 0.;

      return zeta3
         + PI2*zm1/6.
         + (ceil*IPI - lzm1/2. - 1./12.*pow2(3.*I + PI))*pow2(zm1)
         + (lzm1/2. + 1./36.*(-21. - 18.*(-1 + 2*ceil)*IPI + 2*PI2))*pow3(zm1)
         + (-11./24.*lzm1 + 1./288.*(131. + 132.*(-1 + 2*ceil)*IPI - 12*PI2))*pow4(zm1)
         + (5./12.*lzm1 + 1./720.*(-265. - 300.*(-1 + 2*ceil)*IPI + 24*PI2))*pow5(zm1)
         + (-137./360.*lzm1 + 1./7200.*(2213. + 2740.*(-1 + 2*ceil)*IPI - 200*PI2))*pow6(zm1)
         + (-947./3600. - 7./20.*(-1 + 2*ceil)*IPI + 7./20.*lzm1 + PI2/42.)*pow7(zm1)
         + (647707./2822400. + 363.*(-1 + 2*ceil)*IPI/1120. - 363./1120.*lzm1 - PI2/48.)*pow8(zm1);
   }
   if (is_close(z, -1.))
      return -0.75*zeta3;
   if (is_close(z, 0.5)) {
      const double ln2  = 0.6931471805599453; // ln(2)
      const double ln23 = 0.3330246519889295; // ln(2)^3
      return (-2*PI2*ln2 + 4*ln23 + 21*zeta3)/24.;
   }

   std::complex<double> u, sum;

   if (std::abs(z) <= 1.) {
      u = -clog(1. - z);
   } else { // az > 1.
      u = -clog(1. - 1./z);
      sum = -pow3(clog(-z))/6. - PI2/6.*clog(-z);
   }

   std::complex<double> p = 1.;

   for (const double b: bf) {
      p *= u;
      sum += b*p;
   }

   return sum;
}

} // namespace polylogarithm
