/* ====================================================================
 * This file is part of Polylogarithm.
 *
 * Polylogarithm is licenced under the MIT License.
 * ==================================================================== */

#include <complex.h>
#include <float.h>
#include <math.h>
#include "fast_clog.h"


static long double hornerl(long double x, const long double* c, int len)
{
   long double p = 0;
   while (len--)
      p = p*x + c[len];
   return p;
}


/**
 * @brief Real dilogarithm \f$\operatorname{Li}_2(x)\f$
 * @param x real argument
 * @return \f$\operatorname{Li}_2(x)\f$
 * @author Alexander Voigt
 *
 * Implemented as a rational function approximation with a maximum
 * error of 8e-9.
 */
float li2f(float x)
{
   const float PI = 3.14159265f;
   const float P[] = {
      9.99999992e-1f,
     -1.07392585e+0f,
      2.29935934e-1f
   };
   const float Q[] = {
      1.00000000e+0f,
     -1.32392672e+0f,
      4.49822034e-1f,
     -2.79507771e-2f
   };

   float y = 0, r = 0, s = 1;

   /* transform to [0, 1/2] */
   if (x < -1) {
      const float l = logf(1 - x);
      y = 1/(1 - x);
      r = -PI*PI/6 + l*(0.5f*l - logf(-x));
      s = 1;
   } else if (x == -1) {
      return -PI*PI/12;
   } else if (x < 0) {
      const float l = log1pf(-x);
      y = x/(x - 1);
      r = -0.5f*l*l;
      s = -1;
   } else if (x == 0) {
      return 0;
   } else if (x < 0.5f) {
      y = x;
      r = 0;
      s = 1;
   } else if (x < 1) {
      y = 1 - x;
      r = PI*PI/6 - logf(x)*log1pf(-x);
      s = -1;
   } else if (x == 1) {
      return PI*PI/6;
   } else if (x < 2) {
      const float l = logf(x);
      y = 1 - 1/x;
      r = PI*PI/6 - l*(logf(y) + 0.5f*l);
      s = 1;
   } else {
      const float l = logf(x);
      y = 1/x;
      r = PI*PI/3 - 0.5f*l*l;
      s = -1;
   }

   const float y2 = y*y;
   const float p = P[0] + y * P[1] + y2 * P[2];
   const float q = Q[0] + y * Q[1] + y2 * (Q[2] + y * Q[3]);

   return r + s*y*p/q;
}

/**
 * @brief Real dilogarithm \f$\operatorname{Li}_2(x)\f$
 * @param x real argument
 * @return \f$\operatorname{Li}_2(x)\f$
 * @author Alexander Voigt
 *
 * Implemented as a rational function approximation with a maximum
 * error of 5e-17
 * [[arXiv:2201.01678](https://arxiv.org/abs/2201.01678)].
 */
double li2(double x)
{
   const double PI = 3.1415926535897932;
   const double P[] = {
      0.9999999999999999502e+0,
     -2.6883926818565423430e+0,
      2.6477222699473109692e+0,
     -1.1538559607887416355e+0,
      2.0886077795020607837e-1,
     -1.0859777134152463084e-2
   };
   const double Q[] = {
      1.0000000000000000000e+0,
     -2.9383926818565635485e+0,
      3.2712093293018635389e+0,
     -1.7076702173954289421e+0,
      4.1596017228400603836e-1,
     -3.9801343754084482956e-2,
      8.2743668974466659035e-4
   };

   double y = 0, r = 0, s = 1;

   /* transform to [0, 1/2] */
   if (x < -1) {
      const double l = log(1 - x);
      y = 1/(1 - x);
      r = -PI*PI/6 + l*(0.5*l - log(-x));
      s = 1;
   } else if (x == -1) {
      return -PI*PI/12;
   } else if (x < 0) {
      const double l = log1p(-x);
      y = x/(x - 1);
      r = -0.5*l*l;
      s = -1;
   } else if (x == 0) {
      return 0;
   } else if (x < 0.5) {
      y = x;
      r = 0;
      s = 1;
   } else if (x < 1) {
      y = 1 - x;
      r = PI*PI/6 - log(x)*log1p(-x);
      s = -1;
   } else if (x == 1) {
      return PI*PI/6;
   } else if (x < 2) {
      const double l = log(x);
      y = 1 - 1/x;
      r = PI*PI/6 - l*(log(y) + 0.5*l);
      s = 1;
   } else {
      const double l = log(x);
      y = 1/x;
      r = PI*PI/3 - 0.5*l*l;
      s = -1;
   }

   const double y2 = y*y;
   const double y4 = y2*y2;
   const double p = P[0] + y * P[1] + y2 * (P[2] + y * P[3]) +
                    y4 * (P[4] + y * P[5]);
   const double q = Q[0] + y * Q[1] + y2 * (Q[2] + y * Q[3]) +
                    y4 * (Q[4] + y * Q[5] + y2 * Q[6]);

   return r + s*y*p/q;
}


/**
 * @brief Real dilogarithm \f$\operatorname{Li}_2(z)\f$ with long double precision
 * @param x real argument
 * @return \f$\operatorname{Li}_2(z)\f$
 * @author Alexander Voigt
 *
 * Implemented as an economized Pade approximation with a maximum
 * error of 2.13e-20 (long double precision) and 1.03e-38 (quadruple
 * precision), respectively.
 */
long double li2l(long double x)
{
   const long double PI  = 3.14159265358979323846264338327950288L;

#if LDBL_DIG <= 18
   const long double P[] = {
      1.07061055633093042767673531395124630e+0L,
     -5.25056559620492749887983310693176896e+0L,
      1.03934845791141763662532570563508185e+1L,
     -1.06275187429164237285280053453630651e+1L,
      5.95754800847361224707276004888482457e+0L,
     -1.78704147549824083632603474038547305e+0L,
      2.56952343145676978700222949739349644e-1L,
     -1.33237248124034497789318026957526440e-2L,
      7.91217309833196694976662068263629735e-5L
   };
   const long double Q[] = {
      1.00000000000000000000000000000000000e+0L,
     -5.20360694854541370154051736496901638e+0L,
      1.10984640257222420881180161591516337e+1L,
     -1.24997590867514516374467903875677930e+1L,
      7.97919868560471967115958363930214958e+0L,
     -2.87732383715218390800075864637472768e+0L,
      5.49210416881086355164851972523370137e-1L,
     -4.73366369162599860878254400521224717e-2L,
      1.23136575793833628711851523557950417e-3L
   };
#else
   const long double P[] = {
      1.0706105563309304276767353139512463033e+0L,
     -1.0976999936495889694316759019749736514e+1L,
      5.0926847456521671435598196295391041399e+1L,
     -1.4137651316583179225403308056234710093e+2L,
      2.6167438966024905838946321639772104746e+2L,
     -3.4058684964478756354428571563597134889e+2L,
      3.2038320769322335817541819891481680235e+2L,
     -2.2042602382067574460998180123431383994e+2L,
      1.1098617775200190748144566197068215051e+2L,
     -4.0511655788875306581280201485439376229e+1L,
      1.0505482055068989911823839054507173824e+1L,
     -1.8711701614474515075977773430379529915e+0L,
      2.1700009060344649725726270289240541652e-1L,
     -1.5030137964245010798524220531488795883e-2L,
      5.3441650657913964794668981233490297624e-4L,
     -7.0813883289232115572593407126124042503e-6L,
      9.9037749983459843207756017329825486520e-9L
   };
   const long double Q[] = {
      1.0000000000000000000000000000000000000e+0L,
     -1.0552362671556327276432317611336360654e+1L,
      5.0559576537049301822461383847757604795e+1L,
     -1.4554341156550484644439941130107249144e+2L,
      2.8070530313005385434963768448877956414e+2L,
     -3.8295927403204227055447294669165437162e+2L,
      3.8034123327297619837621395664129395552e+2L,
     -2.7878320883603471098333035075019703189e+2L,
      1.5127204204944666467295435798726892255e+2L,
     -6.0401294167088266781532200159750431265e+1L,
      1.7480717870136655299932069333598290551e+1L,
     -3.5731199220918363429317367913920234597e+0L,
      4.9537900060585746351429453828993042241e-1L,
     -4.3750415383481013167382471956794084257e-2L,
      2.2234568413141078158637017098109483650e-3L,
     -5.4149527323164210917629185927908003034e-5L,
      4.1629116823949062782848730828350635143e-7L
   };
#endif

   long double y = 0, r = 0, s = 1;

   /* transform to [0, 1/2) */
   if (x < -1) {
      const long double l = logl(1 - x);
      y = 1/(1 - x);
      r = -PI*PI/6 + l*(0.5L*l - logl(-x));
      s = 1;
   } else if (x == -1) {
      return -PI*PI/12;
   } else if (x < 0) {
      const long double l = log1pl(-x);
      y = x/(x - 1);
      r = -0.5L*l*l;
      s = -1;
   } else if (x == 0) {
      return 0;
   } else if (x < 0.5L) {
      y = x;
      r = 0;
      s = 1;
   } else if (x < 1) {
      y = 1 - x;
      r = PI*PI/6 - logl(x)*log1pl(-x);
      s = -1;
   } else if (x == 1) {
      return PI*PI/6;
   } else if (x < 2) {
      const long double l = logl(x);
      y = 1 - 1/x;
      r = PI*PI/6 - l*(logl(y) + 0.5L*l);
      s = 1;
   } else {
      const long double l = logl(x);
      y = 1/x;
      r = PI*PI/3 - 0.5L*l*l;
      s = -1;
   }

   const long double z = y - 0.25L;

   const long double p = hornerl(z, P, sizeof(P)/sizeof(P[0]));
   const long double q = hornerl(z, Q, sizeof(Q)/sizeof(Q[0]));

   return r + s*y*p/q;
}


/**
 * @brief Complex dilogarithm \f$\operatorname{Li}_2(z)\f$
 * @param z complex argument
 * @return \f$\operatorname{Li}_2(z)\f$
 * @note Implementation translated from SPheno to C++
 * @author Werner Porod
 * @note translated to C++ by Alexander Voigt
 */
double _Complex cli2(double _Complex z)
{
   const double PI = 3.1415926535897932;

   /* bf[1..N-1] are the even Bernoulli numbers / (2 n + 1)! */
   /* generated by: Table[BernoulliB[2 n]/(2 n + 1)!, {n, 1, 9}] */
   const double bf[10] = {
      - 1.0/4.0,
      + 1.0/36.0,
      - 1.0/3600.0,
      + 1.0/211680.0,
      - 1.0/10886400.0,
      + 1.0/526901760.0,
      - 4.0647616451442255e-11,
      + 8.9216910204564526e-13,
      - 1.9939295860721076e-14,
      + 4.5189800296199182e-16
   };

   const double rz = creal(z);
   const double iz = cimag(z);

   /* special cases */
   if (iz == 0.0) {
      if (rz <= 1.0) {
         return li2(rz);
      }
      // rz > 1.0
      return li2(rz) - PI*log(rz)*I;
   }

   const double nz = rz*rz + iz*iz;

   if (nz < DBL_EPSILON) {
      return z*(1.0 + 0.25*z);
   }

   double _Complex u = 0.0, rest = 0.0;
   int sgn = 1;

   /* transformation to |z|<1, Re(z)<=0.5 */
   if (rz <= 0.5) {
      if (nz > 1.0) {
         const double _Complex lz = fast_clog(-z);
         u = -fast_clog(1.0 - 1.0 / z);
         rest = -0.5*lz*lz - PI*PI/6;
         sgn = -1;
      } else { /* nz <= 1 */
         u = -fast_clog(1.0 - z);
         rest = 0;
         sgn = 1;
      }
   } else { /* rz > 0.5 */
      if (nz <= 2*rz) {
         u = -fast_clog(z);
         rest = u*fast_clog(1.0 - z) + PI*PI/6;
         sgn = -1;
      } else { /* nz > 2*rz */
         const double _Complex lz = fast_clog(-z);
         u = -fast_clog(1.0 - 1.0 / z);
         rest = -0.5*lz*lz - PI*PI/6;
         sgn = -1;
      }
   }

   const double _Complex u2 = u*u;
   const double _Complex u4 = u2*u2;
   const double _Complex sum =
      u +
      u2 * (bf[0] +
      u  * (bf[1] +
      u2 * (
          bf[2] +
          u2*bf[3] +
          u4*(bf[4] + u2*bf[5]) +
          u4*u4*(bf[6] + u2*bf[7] + u4*(bf[8] + u2*bf[9]))
      )));

   return (double)sgn * sum + rest;
}


/**
 * @brief Complex dilogarithm \f$\operatorname{Li}_2(z)\f$ with long double precision
 * @param z complex argument
 * @return \f$\operatorname{Li}_2(z)\f$
 * @note Implementation translated from SPheno to C++
 * @author Werner Porod
 * @note translated to C++ and extended to long double precision by Alexander Voigt
 */
long double _Complex cli2l(long double _Complex z)
{
   const long double PI = 3.14159265358979323846264338327950288L;

   /* bf[1..N-1] are the even Bernoulli numbers / (2 n + 1)! */
   /* generated by: Table[BernoulliB[2 n]/(2 n + 1)!, {n, 1, 22}] */
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

   const long double rz = creall(z);
   const long double iz = cimagl(z);

   /* special cases */
   if (iz == 0) {
      if (rz <= 1) {
         return li2l(rz);
      }
      // rz > 1
      return li2l(rz) - PI*logl(rz)*I;
   }

   const long double nz = rz*rz + iz*iz;

   if (nz < LDBL_EPSILON) {
      return z*(1.0L + 0.25L*z);
   }

   long double _Complex u = 0, rest = 0;
   int sgn = 1;

   /* transformation to |z|<1, Re(z)<=0.5 */
   if (rz <= 0.5L) {
      if (nz > 1.0L) {
         const long double _Complex lz = fast_clogl(-z);
         u = -fast_clogl(1.0L - 1.0L/z);
         rest = -0.5L*lz*lz - PI*PI/6;
         sgn = -1;
      } else { /* nz <= 1 */
         u = -fast_clogl(1.0L - z);
         rest = 0;
         sgn = 1;
      }
   } else { /* rz > 0.5L */
      if (nz <= 2*rz) {
         u = -fast_clogl(z);
         rest = u * fast_clogl(1.0L - z) + PI*PI/6;
         sgn = -1;
      } else { /* nz > 2*rz */
         const long double _Complex lz = fast_clogl(-z);
         u = -fast_clogl(1.0L - 1.0L/z);
         rest = -0.5L*lz*lz - PI*PI/6;
         sgn = -1;
      }
   }

   /* the dilogarithm */
   const long double _Complex u2 = u*u;
   long double _Complex sum = 0.0L;

   for (int i = sizeof(bf)/sizeof(bf[0]) - 1; i >= 2; i--) {
      sum = u2 * (bf[i] + sum);
   }

   /* lowest order terms w/ different powers */
   sum = u + u2 * (bf[0] + u * (bf[1] + sum));

   return (long double)sgn * sum + rest;
}
