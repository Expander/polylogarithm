// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the MIT License.
// ====================================================================

#include "Li5.hpp"
#include "complex.hpp"
#include <cfloat>
#include <cmath>

namespace polylogarithm {

namespace {

   template <typename T, int N>
   Complex<T> horner(const Complex<T>& z, const T (&coeffs)[N]) noexcept
   {
      static_assert(N >= 2, "more than two coefficients required");

      const T r = z.re + z.re;
      const T s = z.re * z.re + z.im * z.im;
      T a = coeffs[N - 1], b = coeffs[N - 2];

      for (int i = N - 3; i >= 0; --i) {
         const T t = a;
         a = b + r * a;
         b = coeffs[i] - s * t;
      }

      return Complex<T>(z.re*a + b, z.im*a);
   }

} // anonymous namespace

/**
 * @brief Complex polylogarithm \f$\operatorname{Li}_5(z)\f$
 * @param z_ complex argument
 * @return \f$\operatorname{Li}_5(z)\f$
 * @author Alexander Voigt
 */
std::complex<double> Li5(const std::complex<double>& z_) noexcept
{
   const double PI    = 3.1415926535897932;
   const double PI2   = PI*PI;
   const double PI4   = PI2*PI2;
   const double zeta5 = 1.0369277551433699;
   const double bf[19] = {
      1.0                   , -15.0/32.0             ,
      1.3953189300411523e-01, -2.8633777006172840e-02,
      4.0317412551440329e-03, -3.3985018004115226e-04,
      4.5445184621617666e-06,  2.3916808048569012e-06,
     -1.2762692600122747e-07, -3.1628984306505932e-08,
      3.2848118445335192e-09,  4.7613713995660579e-10,
     -8.0846898171909830e-11, -7.2387648587737207e-12,
      1.9439760115173968e-12,  1.0256978405977236e-13,
     -4.6180551009884830e-14, -1.1535857196470580e-15,
      1.0903545401333394e-15
   };

   const Complex<double> z = { std::real(z_), std::imag(z_) };

   if (z.im == 0) {
      if (z.re == 0) {
         return 0.0;
      }
      if (z.re == 1) {
         return zeta5;
      }
      if (z.re == -1) {
         return -15.0*zeta5/16.0;
      }
   }

   const double nz  = norm_sqr(z);
   const double pz  = arg(z);
   const double lnz = 0.5*std::log(nz);

   if (lnz*lnz + pz*pz < 1) { // |log(z)| < 1
      const Complex<double> u(lnz, pz); // log(z)
      const Complex<double> u2 = u*u;
      const double c0 = zeta5;
      const double c1 = 1.0823232337111382; // zeta(4)
      const double c2 = 0.60102845157979714; // zeta(3)/2
      const double c3 = 0.27415567780803774;
      const Complex<double> c4 = (25.0/12.0 - log(-u))/24.0;
      const double c5 = -1.0/240.0;

      const double cs[6] = {
         -1.1574074074074074e-04, 2.0667989417989418e-07,
         -1.0935444136502338e-09, 8.6986487449450412e-12,
         -8.6899587861588824e-14, 1.0081254080218813e-15
      };

      return c0 + u * c1 +
         u2 * (c2 + u * c3 +
         u2 * (c4 + u * c5 +
         u2 * horner(u2, cs)));
   }

   Complex<double> u(0.0, 0.0), rest(0.0, 0.0);

   if (nz <= 1) {
      u = -log(1.0 - z);
   } else { // nz > 1
      const double arg = pz > 0.0 ? pz - PI : pz + PI;
      const Complex<double> lmz(lnz, arg); // log(-z)
      const Complex<double> lmz2 = lmz*lmz;
      u = -log(1.0 - 1.0/z);
      rest = -1.0/360.0*lmz*(7*PI4 + lmz2*(10.0*PI2 + 3.0*lmz2));
   }

   const Complex<double> u2 = u*u;
   const Complex<double> u4 = u2*u2;
   const Complex<double> u8 = u4*u4;

   return
      rest +
      u*bf[0] +
      u2*(bf[1] + u*bf[2]) +
      u4*(bf[3] + u*bf[4] + u2*(bf[5] + u*bf[6])) +
      u8*(bf[7] + u*bf[8] + u2*(bf[9] + u*bf[10]) +
          u4*(bf[11] + u*bf[12] + u2*(bf[13] + u*bf[14]))) +
      u8*u8*(bf[15] + u*bf[16] + u2*(bf[17] + u*bf[18]));
}

/**
 * @brief Complex polylogarithm \f$\operatorname{Li}_5(z)\f$ with long double precision
 * @param z_ complex argument
 * @return \f$\operatorname{Li}_5(z)\f$
 * @author Alexander Voigt
 */
std::complex<long double> Li5(const std::complex<long double>& z_) noexcept
{
   const long double PI    = 3.14159265358979323846264338327950288L;
   const long double PI2   = PI*PI;
   const long double PI4   = PI2*PI2;
   const long double zeta5 = 1.03692775514336992633136548645703417L;
   const long double bf[] = {
      1.0L,
     -15.0L/32.0L,
      1.39531893004115226337448559670781893e-01L,
     -2.86337770061728395061728395061728395e-02L,
      4.03174125514403292181069958847736626e-03L,
     -3.39850180041152263374485596707818930e-04L,
      4.54451846216176664909446003743133026e-06L,
      2.39168080485690118829088702752453967e-06L,
     -1.27626926001227465885443518981747128e-07L,
     -3.16289843065059324402567872795470007e-08L,
      3.28481184453351916215185742384719818e-09L,
      4.76137139956605790483191328977324967e-10L,
     -8.08468981719098302564602603623317490e-11L,
     -7.23876485877372069468292158375897580e-12L,
      1.94397601151739684930556492093599766e-12L,
      1.02569784059772359718813559433198981e-13L,
     -4.61805510098848301805820862410656158e-14L,
     -1.15358571964705800368425114834190972e-15L,
      1.09035454013333939879770883809662643e-15L,
      2.31481363172925263940797103190091493e-18L,
     -2.56699170432652921943348919933966693e-17L,
#if LDBL_DIG > 18
      4.57086206073149690144959626860139115e-19L,
      6.03667796132057058823561033114107090e-19L,
     -2.16776249440624129587941717218396578e-20L,
     -1.41940966156001652983322668820112130e-20L,
      7.50200095064138625532377521619527234e-22L,
      3.33870453950783971643715159254469304e-22L,
     -2.30600404426203476825215151352586388e-23L,
     -7.85817324568948189044990646315027350e-24L,
      6.66834530437388085486513704613056895e-25L,
      1.85091565409252971894649796883651603e-25L,
     -1.85915294451740855841031840576364891e-26L,
     -4.36297464803458904472660817437095794e-27L,
      5.06110760995292844822634895878109349e-28L,
      1.02919182497568782037888979300008731e-28L,
     -1.35513912210183166156877853283041765e-29L,
     -2.42940596129573826559241540956570701e-30L,
      3.58519739665037052115164738053066333e-31L,
      5.73796581610397206400846538673280837e-32L,
     -9.40035936245687345352774520389480989e-33L,
     -1.35590280493486311500090284171223034e-33L,
      2.44784384191528918377859141748058708e-34L,
      3.20528849130720958124170712217928968e-35L,
     -6.33983878185254827964152375942174520e-36L,
     -7.57925545801218291534870941639851689e-37L
#endif
   };

   const Complex<long double> z = { std::real(z_), std::imag(z_) };

   if (z.im == 0) {
      if (z.re == 0) {
         return 0.0L;
      }
      if (z.re == 1) {
         return zeta5;
      }
      if (z.re == -1) {
         return -15.0L*zeta5/16.0L;
      }
   }

   const long double nz  = norm_sqr(z);
   const long double pz  = arg(z);
   const long double lnz = 0.5L*std::log(nz);

   if (lnz*lnz + pz*pz < 1) { // |log(z)| < 1
      const Complex<long double> u(lnz, pz); // log(z)
      const Complex<long double> u2 = u*u;
      const long double c0 = zeta5;
      const long double c1 = 1.08232323371113819151600369654116790L; // zeta(4)
      const long double c2 = 0.601028451579797142699869080755724995L; // zeta(3)/2
      const long double c3 = 0.274155677808037739412069194441004198L;
      const Complex<long double> c4 = (25.0L/12.0L - log(-u))/24.0L;
      const long double c5 = -1.0L/240.0L;

      const long double cs[] = {
        -1.15740740740740740740740740740740741e-04L,
         2.06679894179894179894179894179894180e-07L,
        -1.09354441365023375605386187396769407e-09L,
         8.69864874494504124133753763383393013e-12L,
        -8.68995878615888235897855907475917096e-14L,
         1.00812540802188133105305292318749173e-15L,
        -1.30160058071551887185136369583811112e-17L,
#if LDBL_DIG > 18
         1.82193858376471817317584461603964703e-19L,
        -2.71703945984843566116551019291461834e-21L,
         4.26404710646461092493552846764602135e-23L,
        -6.97907523609028772068847244464155123e-25L,
         1.18322350137468824961908885022923872e-26L,
        -2.06699310884539801347149709357488995e-28L,
         3.70514089192917181888714449508744051e-30L,
        -6.79216396752655546304766934905316826e-32L,
         1.26987457489641744061409835124861589e-33L,
        -2.41590408788077922327391223699041147e-35L,
         4.66812326074215685597260023169508118e-37L
#endif
      };

      return c0 + u * c1 +
         u2 * (c2 + u * c3 +
         u2 * (c4 + u * c5 +
         u2 * horner(u2, cs)));
   }

   Complex<long double> u(0.0L, 0.0L), rest(0.0L, 0.0L);

   if (nz <= 1) {
      u = -log(1.0L - z);
   } else { // nz > 1
      const long double arg = pz > 0.0 ? pz - PI : pz + PI;
      const Complex<long double> lmz(lnz, arg); // log(-z)
      const Complex<long double> lmz2 = lmz*lmz;
      u = -log(1.0L - 1.0L/z);
      rest = -1.0L/360.0L*lmz*(7*PI4 + lmz2*(10.0L*PI2 + 3.0L*lmz2));
   }

   return rest + u*horner(u, bf);
}

} // namespace polylogarithm
