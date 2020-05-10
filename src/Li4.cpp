// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the GNU Lesser General Public
// License (GNU LGPL) version 3.
// ====================================================================

#include "Li4.hpp"
#include <cfloat>
#include <cmath>

namespace polylogarithm {

namespace {
   template <typename T>
   std::complex<T> clog(std::complex<T> z) noexcept
   {
      // converts -0.0 to 0.0
      const T rz = std::real(z) == T(0) ? std::abs(std::real(z)) : std::real(z);
      const T iz = std::imag(z) == T(0) ? std::abs(std::imag(z)) : std::imag(z);

      return std::complex<T>(0.5*std::log(rz*rz + iz*iz), std::atan2(iz, rz));
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

   const auto rz  = std::real(z);
   const auto iz  = std::imag(z);

   if (iz == 0) {
      if (rz == 0) {
         return 0.0;
      }
      if (rz == 1) {
         return zeta4;
      }
      if (rz == -1.0) {
         return -7.0*PI4/720.0;
      }
   }

   const auto nz  = rz*rz + iz*iz;
   const auto pz  = std::atan2(iz, rz);
   const auto lnz = 0.5*std::log(nz);

   if (lnz*lnz + pz*pz < 1.0) { // |log(z)| < 1
      const auto u  = std::complex<double>(lnz, pz); // clog(z)
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

      return
         cadd(zeta4,
         cadd(cmul(u2, cadd(c2, cmul(u2, c4))),
            cmul(u,
               cadd(c1,
                  cmul(u2, cadd(c3,
                  cmul(u2, cadd(cs[0],
                  cmul(u2, cadd(cs[1],
                  cmul(u2, cadd(cs[2],
                  cmul(u2, cadd(cs[3],
                  cmul(u2, cadd(cs[4],
                  cmul(u2, cadd(cs[5],
                  cmul(u2, cs[6])))))))))))))))
                  )
               )
            )
         );
   }

   std::complex<double> u(0.0, 0.0), r(0.0, 0.0);
   double sgn = 1;

   if (nz <= 1.0) {
      u = -clog(1.0 - z);
   } else { // nz > 1
      const auto arg = pz > 0.0 ? pz - PI : pz + PI;
      const auto lmz = std::complex<double>(lnz, arg); // clog(-z)
      const auto lmz2 = lmz*lmz;
      u = -clog(1.0 - 1.0/z);
      r = 1.0/360.0*(-7*PI4 + lmz2*(-30.0*PI2 - 15.0*lmz2));
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

/**
 * @brief Complex polylogarithm \f$\mathrm{Li}_4(z)\f$ with long double precision
 * @param z complex argument
 * @return \f$\mathrm{Li}_4(z)\f$
 */
std::complex<long double> Li4(const std::complex<long double>& z)
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

   const auto rz  = std::real(z);
   const auto iz  = std::imag(z);

   if (iz == 0) {
      if (rz == 0) {
         return 0.0L;
      }
      if (rz == 1) {
         return zeta4;
      }
      if (rz == -1.0) {
         return -7.0L*PI4/720.0L;
      }
   }

   const auto nz  = rz*rz + iz*iz;
   const auto pz  = std::atan2(iz, rz);
   const auto lnz = 0.5L*std::log(nz);

   if (lnz*lnz + pz*pz < 1.0L) { // |log(z)| < 1
      const auto u  = std::complex<long double>(lnz, pz); // clog(z)
      const auto u2 = u*u;
      const auto c1 = 1.20205690315959428539973816151144999L; // zeta(3)
      const auto c2 = 0.822467033424113218236207583323012595L;
      const auto c3 = (11.0L/6.0L - clog(-u))/6.0L;
      const auto c4 = -1.0L/48.0L;

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

      std::complex<long double> sum(0.0L, 0.0L);

      for (int i = sizeof(cs)/sizeof(cs[0]) - 1; i >= 0; i--) {
         sum = cmul(u2, cadd(cs[i], sum));
      }

      // lowest order terms w/ different powers
      sum = cadd(zeta4, cadd(cmul(u2, cadd(c2, cmul(u2, c4))),
                             cmul(u, cadd(c1, cmul(u2, cadd(c3, sum))))));

      return sum;
   }

   std::complex<long double> u(0.0L, 0.0L), r(0.0L, 0.0L);
   long double sgn = 1;

   if (nz <= 1.0L) {
      u = -clog(1.0L - z);
   } else { // nz > 1
      const auto arg = pz > 0.0 ? pz - PI : pz + PI;
      const auto lmz = std::complex<long double>(lnz, arg); // clog(-z)
      const auto lmz2 = lmz*lmz;
      u = -clog(1.0L - 1.0L/z);
      r = 1.0L/360.0L*(-7*PI4 + lmz2*(-30.0L*PI2 - 15.0L*lmz2));
      sgn = -1;
   }

   std::complex<long double> sum(0.0L, 0.0L);

   for (int i = sizeof(bf)/sizeof(bf[0]) - 1; i >= 0; i--) {
      sum = cmul(u, cadd(bf[i], sum));
   }

   return sgn*sum + r;
}

} // namespace polylogarithm
