// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the GNU Lesser General Public
// License (GNU LGPL) version 3.
// ====================================================================

#include "Li2.hpp"
#include "complex.hpp"
#include <cfloat>
#include <cmath>
#include <limits>

namespace polylogarithm {

namespace {

   template <typename T>
   T horner(T x, const T* c, int len)
   {
      T p = 0;
      while (len--)
         p = p*x + c[len];
      return p;
   }

   template <int Nstart, int Nend, typename T, int N>
   Complex<T> horner(const Complex<T>& z, const T (&coeffs)[N]) noexcept
   {
      static_assert(Nstart <= Nend && Nend < N && Nend >= 1, "invalid array bounds");

      const T r = z.re + z.re;
      const T s = z.re * z.re + z.im * z.im;
      T a = coeffs[Nend], b = coeffs[Nend - 1];

      for (int i = Nend - 2; i >= Nstart; --i) {
         const T t = a;
         a = b + r * a;
         b = coeffs[i] - s * t;
      }

      return Complex<T>(z.re*a + b, z.im*a);
   }

} // anonymous namespace

/**
 * @brief Real dilogarithm \f$\mathrm{Li}_2(x)\f$
 * @param x real argument
 * @return \f$\mathrm{Li}_2(x)\f$
 * @author Alexander Voigt
 *
 * Implemented as a rational function approximation with a
 * maximum relative error of 2.1e-18.
 */
double Li2(double x) noexcept
{
   const double PI = 3.141592653589793;
   const double P[] = {
      1.0000000000000000021e+0,
     -2.9439077825945964913e+0,
      3.2670342301938645239e+0,
     -1.6822678547319245357e+0,
      3.9417801505613831897e-1,
     -3.3584437418278648078e-2,
      3.6579548569932482700e-4
   };
   const double Q[] = {
      1,
     -3.1939077825945954406e+0,
      3.9544000647313133535e+0,
     -2.3784892284012652398e+0,
      7.0904177354688084044e-1,
     -9.3739769867521022404e-2,
      3.8093869120116787591e-3
   };

   if (x == 0)
      return 0;
   if (x == 1)
      return PI*PI/6;
   if (x == -1)
      return -PI*PI/12;

   double y = 0, r = 0, s = 1;

   /* transform to [0, 1/2] */
   if (x < -1) {
      const double l = std::log(1 - x);
      y = 1/(1 - x);
      r = -PI*PI/6 + l*(0.5*l - std::log(-x));
      s = 1;
   } else if (x < 0) {
      const double l = std::log1p(-x);
      y = x/(x - 1);
      r = -0.5*l*l;
      s = -1;
   } else if (x < 0.5) {
      y = x;
      r = 0;
      s = 1;
   } else if (x < 1) {
      y = 1 - x;
      r = PI*PI/6 - std::log(x)*std::log(1 - x);
      s = -1;
   } else if (x < 2) {
      const double l = std::log(x);
      y = 1 - 1/x;
      r = PI*PI/6 - l*(std::log(1 - 1/x) + 0.5*l);
      s = 1;
   } else {
      const double l = std::log(x);
      y = 1/x;
      r = PI*PI/3 - 0.5*l*l;
      s = -1;
   }

   const double p = horner(y, P, sizeof(P)/sizeof(P[0]));
   const double q = horner(y, Q, sizeof(Q)/sizeof(Q[0]));

   return r + s*y*p/q;
}

/**
 * @brief Real dilogarithm \f$\mathrm{Li}_2(z)\f$ with long double precision
 * @param x real argument
 * @return \f$\mathrm{Li}_2(z)\f$
 * @author Alexander Voigt
 *
 * Implemented as a rational function approximation with a maximum
 * relative error of 4.1e-21 (long double) and 3.1 (quadruple
 * precision).
 */
long double Li2(long double x) noexcept
{
   const long double PI = 3.14159265358979323846264338327950288L;
#if LDBL_DIG <= 18
   const long double P[] = {
      1.000000000000000000004139242230383319580e+0L,
     -3.487762980603640658103611020138020359177e+0L,
      4.806535516830445654339162187560030671816e+0L,
     -3.308613600271641896478346549030223509476e+0L,
      1.179396551716581370152377778000509525433e+0L,
     -2.018700678934015014829496389063651400501e-1L,
      1.287456041852283933114884438200724023392e-2L,
     -1.033218746464088069519331444271701827352e-4L
   };
   const long double Q[] = {
      1.000000000000000000000000000000000000000e+0L,
     -3.737762980603640655360823697458395397077e+0L,
      5.629865150870244400485990049358005325948e+0L,
     -4.363272890144340400922968436762690521699e+0L,
      1.838284388221160440266981075256880968845e+0L,
     -4.067668965223713893864200244855877582547e-1L,
      4.124099939142561145287545934805943708748e-2L,
     -1.308374236546985200155320229304373514715e-3L
   };
#else
   const long double P[] = {
      1.000000000000000000000000000000000000000e+0L,
     -7.287836026293007036950639387021726765663e+0L,
      23.81545623790376778894331890102886403563e+0L,
     -46.05033532842722132118308962006325124290e+0L,
      58.54881768373102592597834455330905597413e+0L,
     -51.43421107893188978781352321887201670161e+0L,
      31.92595848973530035544399399745039290324e+0L,
     -14.06988451302420353649035303317666423752e+0L,
      4.359706642554451028170619355055660927151e+0L,
     -9.261359685568062652651023592559153471944e-1L,
      1.288432888720493539643303789043503094981e-1L,
     -1.085087084358217119158266823617666800942e-2L,
      4.796175136914686220858905279884620056548e-4L,
     -8.170779455874784395134674021307361681514e-6L,
      1.576137369111657930232047018998161178865e-8L
   };
   const long double Q[] = {
      1.000000000000000000000000000000000000000e+0L,
     -7.537836026293007036950639387021726764886e+0L,
      25.58880413336590843706986763667318429579e+0L,
     -51.67249902551391987078937437511796918107e+0L,
      69.05485673248994011825511781580902714637e+0L,
     -64.28210107650314382757288000492042122129e+0L,
      42.71867709425390958221672470482605169080e+0L,
     -20.42871842270475507140771110978257576181e+0L,
      6.994363330607140789641413570234544197271e+0L,
     -1.683917502924445721756542581983598343874e+0L,
      2.756356619926038754088086416897457355184e-1L,
     -2.899348828771590606786384940125401751121e-2L,
      1.781707396338762194616290891089680894814e-3L,
     -5.363089470000456715142977989844950980813e-5L,
      5.247190504505569429431039126190470985163e-7L
   };
#endif

   if (x == 0)
      return 0;
   if (x == 1)
      return PI*PI/6;
   if (x == -1)
      return -PI*PI/12;

   long double y = 0, r = 0, s = 1;

   /* transform to [0, 1/2] */
   if (x < -1) {
      const long double l = std::log(1 - x);
      y = 1/(1 - x);
      r = -PI*PI/6 + l*(0.5L*l - std::log(-x));
      s = 1;
   } else if (x < 0) {
      const long double l = std::log1p(-x);
      y = x/(x - 1);
      r = -0.5L*l*l;
      s = -1;
   } else if (x < 0.5L) {
      y = x;
      r = 0;
      s = 1;
   } else if (x < 1) {
      y = 1 - x;
      r = PI*PI/6 - std::log(x)*std::log(1 - x);
      s = -1;
   } else if (x < 2) {
      const long double l = std::log(x);
      y = 1 - 1/x;
      r = PI*PI/6 - l*(std::log(1 - 1/x) + 0.5L*l);
      s = 1;
   } else {
      const long double l = std::log(x);
      y = 1/x;
      r = PI*PI/3 - 0.5L*l*l;
      s = -1;
   }

   const long double p = horner(y, P, sizeof(P)/sizeof(P[0]));
   const long double q = horner(y, Q, sizeof(Q)/sizeof(Q[0]));

   return r + s*y*p/q;
}

/**
 * @brief Complex dilogarithm \f$\mathrm{Li}_2(z)\f$
 * @param z_ complex argument
 * @return \f$\mathrm{Li}_2(z)\f$
 * @note Implementation translated from SPheno to C++
 * @author Werner Porod
 * @note translated to C++ by Alexander Voigt
 */
std::complex<double> Li2(const std::complex<double>& z_) noexcept
{
   const double PI = 3.141592653589793;
   const Complex<double> z = { std::real(z_), std::imag(z_) };

   // bf[1..N-1] are the even Bernoulli numbers / (2 n + 1)!
   // generated by: Table[BernoulliB[2 n]/(2 n + 1)!, {n, 1, 9}]
   const double bf[10] = {
      - 1.0/4.0,
      + 1.0/36.0,
      - 1.0/3600.0,
      + 1.0/211680.0,
      - 1.0/10886400.0,
      + 1.0/526901760.0,
      - 4.064761645144226e-11,
      + 8.921691020456453e-13,
      - 1.993929586072108e-14,
      + 4.518980029619918e-16
   };

   const double nz = norm_sqr(z);

   // special cases
   if (z.im == 0) {
      if (z.re <= 1) {
         return Li2(z.re);
      }
      // z.re > 1
      return { Li2(z.re), -PI*std::log(z.re) };
   } else if (nz < std::numeric_limits<double>::epsilon()) {
      return z_;
   }

   Complex<double> cy(0.0, 0.0), cz(0.0, 0.0);
   double sgn = 1;

   // transformation to |z|<1, Re(z)<=0.5
   if (z.re <= 0.5) {
      if (nz > 1) {
         const Complex<double> lz = log(-z);
         cy = -0.5 * lz*lz - PI * PI / 6.0;
         cz = -log(1.0 - 1.0 / z);
         sgn = -1;
      } else { // nz <= 1
         cy = 0;
         cz = -log(1.0 - z);
         sgn = 1;
      }
   } else { // z.re > 0.5
      if (nz <= 2*z.re) {
         cz = -log(z);
         cy = cz * log(1.0 - z) + PI * PI / 6.0;
         sgn = -1;
      } else { // nz > 2*z.re
         const Complex<double> lz = log(-z);
         cy = -0.5 * lz*lz - PI * PI / 6.0;
         cz = -log(1.0 - 1.0 / z);
         sgn = -1;
      }
   }

   const Complex<double> cz2(cz*cz);

   return sgn*(cz + cz2*(bf[0] + cz*horner<1, 9>(cz2, bf))) + cy;
}

/**
 * @brief Complex dilogarithm \f$\mathrm{Li}_2(z)\f$ with long double precision
 * @param z_ complex argument
 * @return \f$\mathrm{Li}_2(z)\f$
 * @note Implementation translated from SPheno to C++
 * @author Werner Porod
 * @note translated to C++ and extended to long double precision by Alexander Voigt
 */
std::complex<long double> Li2(const std::complex<long double>& z_) noexcept
{
   const long double PI = 3.14159265358979323846264338327950288L;
   const Complex<long double> z = { std::real(z_), std::imag(z_) };

   // bf[1..N-1] are the even Bernoulli numbers / (2 n + 1)!
   // generated by: Table[BernoulliB[2 n]/(2 n + 1)!, {n, 1, 22}]
   const long double bf[] = {
      -1.0L/4.0L                                 ,
       1.0L/36.0L                                ,
      -1.0L/3600.0L                              ,
       1.0L/211680.0L                            ,
      -1.0L/10886400.0L                          ,
       1.0L/526901760.0L                         ,
      -4.06476164514422552680590938629196667e-11L,
       8.92169102045645255521798731675274885e-13L,
      -1.99392958607210756872364434779378971e-14L,
       4.51898002961991819165047655285559323e-16L,
      -1.03565176121812470144834115422186567e-17L,
#if LDBL_DIG > 18
       2.39521862102618674574028374300098038e-19L,
      -5.58178587432500933628307450562541991e-21L,
       1.30915075541832128581230739918659230e-22L,
      -3.08741980242674029324227976486646243e-24L,
       7.31597565270220342035790560925214859e-26L,
      -1.74084565723400074098905514775970255e-27L,
       4.15763564461389971961789962077522667e-29L,
      -9.96214848828462210319400670245583885e-31L,
       2.39403442489616530052116798789374956e-32L,
      -5.76834735536739008429179316187765424e-34L,
       1.39317947964700797782788660391154833e-35L,
      -3.37212196548508947046847363525493096e-37L
#endif
   };

   constexpr int N = sizeof(bf)/sizeof(bf[0]);

   const long double nz = norm_sqr(z);

   // special cases
   if (z.im == 0) {
      if (z.re <= 1) {
         return Li2(z.re);
      }
      if (z.re > 1) {
         return { Li2(z.re), -PI*std::log(z.re) };
      }
   } else if (nz < std::numeric_limits<long double>::epsilon()) {
      return z_;
   }

   Complex<long double> cy(0.0L, 0.0L), cz(0.0L, 0.0L);
   long double sgn = 1;

   // transformation to |z|<1, Re(z)<=0.5
   if (z.re <= 0.5L) {
      if (nz > 1) {
         const Complex<long double> lz = log(-z);
         cy = -0.5L * lz*lz - PI * PI / 6.0L;
         cz = -log(1.0L - 1.0L/z);
         sgn = -1;
      } else { // nz <= 1
         cy = 0;
         cz = -log(1.0L - z);
         sgn = 1;
      }
   } else { // z.re > 0.5L
      if (nz <= 2*z.re) {
         cz = -log(z);
         cy = cz * log(1.0L - z) + PI * PI / 6.0L;
         sgn = -1;
      } else { // nz > 2*z.re
         const Complex<long double> lz = log(-z);
         cy = -0.5L * lz*lz - PI * PI / 6.0L;
         cz = -log(1.0L - 1.0L/z);
         sgn = -1;
      }
   }

   const Complex<long double> cz2(cz*cz);

   return sgn*(cz + cz2*(bf[0] + cz*horner<1, N-1>(cz2, bf))) + cy;
}

} // namespace polylogarithm
