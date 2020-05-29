// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the GNU Lesser General Public
// License (GNU LGPL) version 3.
// ====================================================================

#include "Li4.hpp"
#include "complex.hpp"
#include <cfloat>
#include <cmath>

namespace polylogarithm {

namespace {

   template <typename T, int N>
   Complex<T> horner(const Complex<T>& z, const T (&coeffs)[N]) noexcept
   {
      static_assert(N >= 2, "more than two coefficients required");

      const T x = z.re;
      const T y = z.im;
      const T r = x + x;
      const T s = x * x + y * y;
      T a = coeffs[N - 1], b = coeffs[N - 2];

      for (int i = N - 3; i >= 0; --i) {
         const T t = a;
         a = b + r * a;
         b = coeffs[i] - s * t;
      }

      return Complex<T>(x*a + b, y*a);
   }

} // anonymous namespace

/**
 * @brief Complex polylogarithm \f$\mathrm{Li}_4(z)\f$
 * @param z_ complex argument
 * @return \f$\mathrm{Li}_4(z)\f$
 * @author Alexander Voigt
 */
std::complex<double> Li4(const std::complex<double>& z_) noexcept
{
   const double PI    = 3.141592653589793;
   const double PI2   = PI*PI;
   const double PI4   = PI2*PI2;
   const double zeta4 = 1.082323233711138;
   const double bf[18] = {
      1.0, -7.0/16.0,
      1.165123456790123e-01, -1.982060185185185e-02,
      1.927932098765432e-03, -3.105709876543209e-05,
     -1.562400911485783e-05,  8.485123546773206e-07,
      2.290961660318971e-07, -2.183261421852691e-08,
     -3.882824879172015e-09,  5.446292103220332e-10,
      6.960805210682725e-11, -1.337573768644521e-11,
     -1.278485268526657e-12,  3.260562858024892e-13,
      2.364757116861825e-14, -7.923135122031161e-15,
   };

   const Complex<double> z = { std::real(z_), std::imag(z_) };

   if (z.im == 0) {
      if (z.re == 0) {
         return 0.0;
      }
      if (z.re == 1) {
         return zeta4;
      }
      if (z.re == -1) {
         return -7.0*PI4/720.0;
      }
   }

   const double nz  = norm_sqr(z);
   const double pz  = arg(z);
   const double lnz = 0.5*std::log(nz);

   if (lnz*lnz + pz*pz < 1) { // |log(z)| < 1
      const Complex<double> u(lnz, pz); // log(z)
      const Complex<double> u2 = u*u;
      const Complex<double> u4 = u2*u2;
      const Complex<double> u8 = u4*u4;
      const double c1 = 1.202056903159594; // zeta(3)
      const double c2 = 0.8224670334241132;
      const Complex<double> c3 = (11.0/6.0 - log(-u))/6.0;
      const double c4 = -1.0/48.0;

      const double cs[7] = {
         -6.944444444444444e-04, 1.653439153439153e-06,
         -1.093544413650234e-08, 1.043837849393405e-10,
         -1.216594230062244e-12, 1.613000652835010e-14,
         -2.342881045287934e-16
      };

      return zeta4 + u2 * (c2 + u2 * c4) +
         u * (
             c1 +
             c3*u2 +
             u4*(cs[0] + u2*cs[1]) +
             u8*(cs[2] + u2*cs[3] + u4*(cs[4] + u2*cs[5])) +
             u8*u8*cs[6]
         );
   }

   Complex<double> u(0.0, 0.0), rest(0.0, 0.0);
   double sgn = 1;

   if (nz <= 1) {
      u = -log(1.0 - z);
   } else { // nz > 1
      const double arg = pz > 0.0 ? pz - PI : pz + PI;
      const Complex<double> lmz(lnz, arg); // log(-z)
      const Complex<double> lmz2 = lmz*lmz;
      u = -log(1.0 - 1.0/z);
      rest = 1.0/360.0*(-7*PI4 + lmz2*(-30.0*PI2 - 15.0*lmz2));
      sgn = -1;
   }

   const Complex<double> u2 = u*u;
   const Complex<double> u4 = u2*u2;
   const Complex<double> u8 = u4*u4;

   return
      rest + sgn * (
         u*bf[0] +
         u2*(bf[1] + u*bf[2]) +
         u4*(bf[3] + u*bf[4] + u2*(bf[5] + u*bf[6])) +
         u8*(bf[7] + u*bf[8] + u2*(bf[9] + u*bf[10]) +
             u4*(bf[11] + u*bf[12] + u2*(bf[13] + u*bf[14]))) +
         u8*u8*(bf[15] + u*bf[16] + u2*bf[17])
      );
}

/**
 * @brief Complex polylogarithm \f$\mathrm{Li}_4(z)\f$ with long double precision
 * @param z_ complex argument
 * @return \f$\mathrm{Li}_4(z)\f$
 * @author Alexander Voigt
 */
std::complex<long double> Li4(const std::complex<long double>& z_) noexcept
{
   const long double PI    = 3.14159265358979323846264338327950288L;
   const long double PI2   = PI*PI;
   const long double PI4   = PI2*PI2;
   const long double zeta4 = 1.08232323371113819151600369654116790L;
   const long double bf[] = {
      1.0L,
     -7.0L/16.0L,
      1.16512345679012345679012345679012346e-01L,
     -1.98206018518518518518518518518518519e-02L,
      1.92793209876543209876543209876543210e-03L,
     -3.10570987654320987654320987654320988e-05L,
     -1.56240091148578352983924736435264456e-05L,
      8.48512354677320663715221538350790051e-07L,
      2.29096166031897114453593835470042743e-07L,
     -2.18326142185269169396153523137650122e-08L,
     -3.88282487917201557228066203807765146e-09L,
      5.44629210322033211825798588082320063e-10L,
      6.96080521068272540787723341341208120e-11L,
     -1.33757376864452151995780722036345205e-11L,
     -1.27848526852665716041462463615741700e-12L,
      3.26056285802489224287884181782170918e-13L,
      2.36475711686182573623095048124390137e-14L,
     -7.92313512203116170242999007113724954e-15L,
     -4.34529157099841872504973716264753844e-16L,
      1.92362700625359201161268755267526042e-16L,
      7.81241433319595467072229389687370732e-18L,
     -4.67180384480365552031762824287222012e-18L,
#if LDBL_DIG > 18
     -1.34353443298128478562602226758937243e-19L,
      1.13568268513473432447646983759384846e-19L,
      2.11527562024325868475059834141917946e-21L,
     -2.76420263347465173882817292537310280e-21L,
     -2.70681766082400642561090595581950049e-23L,
      6.73720448286285721432671612656264303e-23L,
      1.32872654566838229758180090450125398e-25L,
     -1.64437730563678264678167631148886630e-24L,
      8.28360589993393411098296734003488096e-27L,
      4.01908484950693506997093150076214959e-26L,
     -4.57571384448487903823597343465369976e-28L,
     -9.83641090946151277583209749821167124e-28L,
      1.69003395560378510677295231219028521e-29L,
      2.41048055630598085046649041649017179e-29L,
     -5.42661270567141825013250340589290005e-31L,
     -5.91424295887417678643375999669283147e-31L,
      1.62321109010873707727111761439681785e-32L,
      1.45275954377402759461325873161579478e-32L,
     -4.65389937002573704417216829815072974e-34L,
     -3.57238626244413318154616242379067282e-34L,
      1.29761714880310295825962542732877943e-35L,
      8.79357407773938851103685229710271214e-36L,
     -3.54800202048240308911663975982519909e-37L
#endif
   };

   const Complex<long double> z = { std::real(z_), std::imag(z_) };

   if (z.im == 0) {
      if (z.re == 0) {
         return 0.0L;
      }
      if (z.re == 1) {
         return zeta4;
      }
      if (z.re == -1) {
         return -7.0L*PI4/720.0L;
      }
   }

   const long double nz  = norm_sqr(z);
   const long double pz  = arg(z);
   const long double lnz = 0.5L*std::log(nz);

   if (lnz*lnz + pz*pz < 1) { // |log(z)| < 1
      const Complex<long double> u(lnz, pz); // log(z)
      const Complex<long double> u2 = u*u;
      const long double c1 = 1.20205690315959428539973816151144999L; // zeta(3)
      const long double c2 = 0.822467033424113218236207583323012595L;
      const Complex<long double> c3 = (11.0L/6.0L - log(-u))/6.0L;
      const long double c4 = -1.0L/48.0L;

      const long double cs[] = {
        -6.94444444444444444444444444444444444e-04L,
         1.65343915343915343915343915343915344e-06L,
        -1.09354441365023375605386187396769407e-08L,
         1.04383784939340494896050451606007162e-10L,
        -1.21659423006224353025699827046628393e-12L,
         1.61300065283501012968488467709998678e-14L,
        -2.34288104528793396933245465250860001e-16L,
         3.64387716752943634635168923207929405e-18L,
#if LDBL_DIG > 18
        -5.97748681166655845456412242441216036e-20L,
         1.02337130555150662198452683223504513e-21L,
        -1.81455956138347480737900283560680332e-23L,
         3.31302580384912709893344878064186841e-25L,
        -6.20097932653619404041449128072466986e-27L,
         1.18564508541733498204388623842798096e-28L,
        -2.30933574895902885743620757867807721e-30L,
         4.57154846962710278621075406449501719e-32L,
        -9.18043553394696104844086650056356358e-34L,
         1.86724930429686274238904009267803247e-35L,
        -3.84151865355610606630482668147264051e-37L
#endif
      };

      return zeta4 + u2*(c2 + u2*c4) +
         u*(c1 + u2*(c3 + u2*horner(u2, cs)));
   }

   Complex<long double> u(0.0L, 0.0L), rest(0.0L, 0.0L);
   long double sgn = 1;

   if (nz <= 1) {
      u = -log(1.0L - z);
   } else { // nz > 1
      const long double arg = pz > 0.0 ? pz - PI : pz + PI;
      const Complex<long double> lmz(lnz, arg); // log(-z)
      const Complex<long double> lmz2 = lmz*lmz;
      u = -log(1.0L - 1.0L/z);
      rest = 1.0L/360.0L*(-7*PI4 + lmz2*(-30.0L*PI2 - 15.0L*lmz2));
      sgn = -1;
   }

   return rest + sgn*u*horner(u, bf);
}

} // namespace polylogarithm
