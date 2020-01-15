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
   const double eps_d = 10.0*std::numeric_limits<double>::epsilon();
   const long double eps_ld = 10.0L*std::numeric_limits<long double>::epsilon();

   template <typename T> T pow2(T x) noexcept { return x*x; }

   // converts -0.0 to 0.0
   template <typename T>
   std::complex<T> clog(std::complex<T> z) noexcept {
      if (std::real(z) == T(0)) { z.real(T(0)); }
      if (std::imag(z) == T(0)) { z.imag(T(0)); }
      return std::log(z);
   }

   template <typename T>
   bool is_close(const std::complex<T>& a, T b, T eps)
   {
      return std::abs(std::real(a) - b) < eps && std::abs(std::imag(a)) < eps;
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
 * @brief Clausen function \f$\mathrm{Cl}_5(\theta) = \mathrm{Re}(\mathrm{Li}_5(e^{i\theta}))\f$
 * @param x real angle
 * @return \f$\mathrm{Cl}_5(\theta)\f$
 */
double Cl5(double x)
{
   const double PI = 3.141592653589793;
   const std::complex<double> i(0.0, 1.0);

   while (x >= 2*PI) {
      x -= 2*PI;
   }

   while (x < 0.0) {
      x += 2*PI;
   }

   return std::real(Li5(std::exp(i*x)));
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
   const double bf[19] = {
      1.0, -15.0/32.0,
      1.395318930041152e-01, -2.863377700617283e-02,
      4.031741255144032e-03, -3.398501800411522e-04,
      4.544518462161766e-06,  2.391680804856901e-06,
     -1.276269260012274e-07, -3.162898430650593e-08,
      3.284811844533519e-09,  4.761371399566057e-10,
     -8.084689817190983e-11, -7.238764858773720e-12,
      1.943976011517396e-12,  1.025697840597723e-13,
     -4.618055100988483e-14, -1.153585719647058e-15,
      1.090354540133339e-15
   };

   if (is_close(z, 0.0, eps_d)) {
      return { 0.0, 0.0 };
   }
   if (is_close(z, 1.0, eps_d)) {
      return { zeta5, 0.0 };
   }
   if (is_close(z, -1.0, eps_d)) {
      return { -15.0*zeta5/16.0, 0.0 };
   }

   const auto az  = std::abs(z);
   const auto pz  = std::arg(z);
   const auto lnz = std::log(az);

   if (pow2(lnz) + pow2(pz) < 1.0) { // |log(z)| < 1
      const auto u  = std::complex<double>(lnz, pz); // clog(z)
      const auto u2 = u*u;
      const auto c0 = zeta5;
      const auto c1 = 1.082323233711138; // zeta(4)
      const auto c2 = 0.6010284515797971; // zeta(3)/2
      const auto c3 = 0.2741556778080377;
      const auto c4 = (25.0/12.0 - clog(-u))/24.0;
      const auto c5 = -1.0/240.0;

      const double cs[6] = {
         -1.157407407407407e-04, 2.066798941798942e-07,
         -1.093544413650234e-09, 8.698648744945041e-12,
         -8.689958786158882e-14, 1.008125408021881e-15
      };

      return
         cadd(c0,
         cadd(cmul(u, c1),
         cmul(u2, cadd(c2,
         cadd(cmul(u, c3),
         cmul(u2, cadd(c4,
         cadd(cmul(u, c5),
         cmul(u2, cadd(cs[0],
         cmul(u2, cadd(cs[1],
         cmul(u2, cadd(cs[2],
         cmul(u2, cadd(cs[3],
         cmul(u2, cadd(cs[4],
         cmul(u2, cs[5])))))))))))))))))));
   }

   std::complex<double> u(0.0, 0.0), rest(0.0, 0.0);

   if (az <= 1.0) {
      u = -clog(1.0 - z);
   } else { // az > 1
      auto arg = PI + pz;
      if (arg > PI) { arg -= 2*PI; }
      const auto lmz = std::complex<double>(lnz, arg); // clog(-z)
      const auto lmz2 = pow2(lmz);
      u = -clog(1.0 - 1.0/z);
      rest = -1.0/360.0*lmz*(7*PI4 + lmz2*(10.0*PI2 + 3.0*lmz2));
   }

   return
      cadd(rest,
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
      cmul(u, cadd(bf[17],
      cmul(u, bf[18]))))))))))))))))))))))))))))))))))))));
}

/**
 * @brief Complex polylogarithm \f$\mathrm{Li}_5(z)\f$ with long double precision
 * @param z complex argument
 * @return \f$\mathrm{Li}_5(z)\f$
 */
std::complex<long double> Li5(const std::complex<long double>& z)
{
   const long double PI    = 3.14159265358979323846264338327950288L;
   const long double PI2   = PI*PI;
   const long double PI4   = PI2*PI2;
   const long double zeta5 = 1.03692775514336992633136548645703417L;
   const long double bf[45] = {
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
   };

   if (is_close(z, 0.0L, eps_ld)) {
      return { 0.0L, 0.0L };
   }
   if (is_close(z, 1.0L, eps_ld)) {
      return { zeta5, 0.0L };
   }
   if (is_close(z, -1.0L, eps_ld)) {
      return { -15.0L*zeta5/16.0L, 0.0L };
   }

   const auto az  = std::abs(z);
   const auto pz  = std::arg(z);
   const auto lnz = std::log(az);

   if (pow2(lnz) + pow2(pz) < 1.0L) { // |log(z)| < 1
      const auto u  = std::complex<long double>(lnz, pz); // clog(z)
      const auto u2 = u*u;
      const auto c0 = zeta5;
      const auto c1 = 1.08232323371113819151600369654116790L; // zeta(4)
      const auto c2 = 0.601028451579797142699869080755724995L; // zeta(3)/2
      const auto c3 = 0.274155677808037739412069194441004198L;
      const auto c4 = (25.0L/12.0L - clog(-u))/24.0L;
      const auto c5 = -1.0L/240.0L;

      const long double cs[18] = {
        -1.15740740740740740740740740740740741e-04L,
         2.06679894179894179894179894179894180e-07L,
        -1.09354441365023375605386187396769407e-09L,
         8.69864874494504124133753763383393013e-12L,
        -8.68995878615888235897855907475917096e-14L,
         1.00812540802188133105305292318749173e-15L,
        -1.30160058071551887185136369583811112e-17L,
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
      };

      return
         cadd(c0,
         cadd(cmul(u, c1),
         cmul(u2, cadd(c2,
         cadd(cmul(u, c3),
         cmul(u2, cadd(c4,
         cadd(cmul(u, c5),
         cmul(u2, cadd(cs[0],
         cmul(u2, cadd(cs[1],
         cmul(u2, cadd(cs[2],
         cmul(u2, cadd(cs[3],
         cmul(u2, cadd(cs[4],
         cmul(u2, cadd(cs[5],
         cmul(u2, cadd(cs[6],
         cmul(u2, cadd(cs[7],
         cmul(u2, cadd(cs[8],
         cmul(u2, cadd(cs[9],
         cmul(u2, cadd(cs[10],
         cmul(u2, cadd(cs[11],
         cmul(u2, cadd(cs[12],
         cmul(u2, cadd(cs[13],
         cmul(u2, cadd(cs[14],
         cmul(u2, cadd(cs[15],
         cmul(u2, cadd(cs[16],
         cmul(u2, cs[17])))))))))))))))))))))))))))))))))))))))))));
   }

   std::complex<long double> u(0.0L, 0.0L), rest(0.0L, 0.0L);

   if (az <= 1.0L) {
      u = -clog(1.0L - z);
   } else { // az > 1
      auto arg = PI + pz;
      if (arg > PI) { arg -= 2*PI; }
      const auto lmz = std::complex<long double>(lnz, arg); // clog(-z)
      const auto lmz2 = pow2(lmz);
      u = -clog(1.0L - 1.0L/z);
      rest = -1.0L/360.0L*lmz*(7*PI4 + lmz2*(10.0L*PI2 + 3.0L*lmz2));
   }

   return
      cadd(rest,
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
      cmul(u, cadd(bf[17],
      cmul(u, cadd(bf[18],
      cmul(u, cadd(bf[19],
      cmul(u, cadd(bf[20],
      cmul(u, cadd(bf[21],
      cmul(u, cadd(bf[22],
      cmul(u, cadd(bf[23],
      cmul(u, cadd(bf[24],
      cmul(u, cadd(bf[25],
      cmul(u, cadd(bf[26],
      cmul(u, cadd(bf[27],
      cmul(u, cadd(bf[28],
      cmul(u, cadd(bf[29],
      cmul(u, cadd(bf[30],
      cmul(u, cadd(bf[31],
      cmul(u, cadd(bf[32],
      cmul(u, cadd(bf[33],
      cmul(u, cadd(bf[34],
      cmul(u, cadd(bf[35],
      cmul(u, cadd(bf[36],
      cmul(u, cadd(bf[37],
      cmul(u, cadd(bf[38],
      cmul(u, cadd(bf[39],
      cmul(u, cadd(bf[40],
      cmul(u, cadd(bf[41],
      cmul(u, cadd(bf[42],
      cmul(u, cadd(bf[43],
      cmul(u, bf[44]))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))));
}

} // namespace polylogarithm
