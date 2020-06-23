/* ====================================================================
 * This file is part of Polylogarithm.
 *
 * Polylogarithm is licenced under the GNU Lesser General Public
 * License (GNU LGPL) version 3.
 * ==================================================================== */

#include <complex.h>
#include <float.h>
#include <math.h>


static double _Complex fast_clog(double _Complex z)
{
   const double rz = creal(z);
   const double iz = cimag(z);

   return 0.5*log(rz*rz + iz*iz) + I*atan2(iz, rz);
}


static long double _Complex fast_clogl(long double _Complex z)
{
   const long double rz = creall(z);
   const long double iz = cimagl(z);

   return 0.5*logl(rz*rz + iz*iz) + I*atan2l(iz, rz);
}


static double horner(double x, const double* c, int len)
{
   double p = 0;
   while (len--)
      p = p*x + c[len];
   return p;
}


static long double hornerl(long double x, const long double* c, int len)
{
   long double p = 0;
   while (len--)
      p = p*x + c[len];
   return p;
}


/**
 * @brief Real dilogarithm \f$\mathrm{Li}_2(x)\f$
 * @param x real argument
 * @return \f$\mathrm{Li}_2(x)\f$
 * @author Alexander Voigt
 *
 * Implemented as an economized Pade approximation with a
 * maximum error of 4.16e-18.
 */
double li2(double x)
{
   const double PI = 3.141592653589793;
   const double P[] = {
      1.0706105563309304277e+0,
     -4.5353562730201404017e+0,
      7.4819657596286408905e+0,
     -6.0516124315132409155e+0,
      2.4733515209909815443e+0,
     -4.6937565143754629578e-1,
      3.1608910440687221695e-2,
     -2.4630612614645039828e-4
   };
   const double Q[] = {
      1.0000000000000000000e+0,
     -4.5355682121856044935e+0,
      8.1790029773247428573e+0,
     -7.4634190853767468810e+0,
      3.6245392503925290187e+0,
     -8.9936784740041174897e-1,
      9.8554565816757007266e-2,
     -3.2116618742475189569e-3
   };

   double y = 0, r = 0, s = 1;

   /* transform to [0, 1/2) */
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
      r = PI*PI/6 - log(x)*log(1 - x);
      s = -1;
   } else if (x == 1) {
      return PI*PI/6;
   } else if (x < 2) {
      const double l = log(x);
      y = 1 - 1/x;
      r = PI*PI/6 - l*(log(1 - 1/x) + 0.5*l);
      s = 1;
   } else {
      const double l = log(x);
      y = 1/x;
      r = PI*PI/3 - 0.5*l*l;
      s = -1;
   }

   const double z = y - 0.25;

   const double p = horner(z, P, sizeof(P)/sizeof(P[0]));
   const double q = horner(z, Q, sizeof(Q)/sizeof(Q[0]));

   return r + s*y*p/q;
}


/**
 * @brief Real dilogarithm \f$\mathrm{Li}_2(z)\f$ with long double precision
 * @param x real argument
 * @return \f$\mathrm{Li}_2(z)\f$
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
      r = PI*PI/6 - logl(x)*logl(1 - x);
      s = -1;
   } else if (x == 1) {
      return PI*PI/6;
   } else if (x < 2) {
      const long double l = logl(x);
      y = 1 - 1/x;
      r = PI*PI/6 - l*(logl(1 - 1/x) + 0.5L*l);
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
 * @brief Complex dilogarithm \f$\mathrm{Li}_2(z)\f$
 * @param z complex argument
 * @note Implementation translated from SPheno to C++
 * @return \f$\mathrm{Li}_2(z)\f$
 */
static double _Complex cli2(double _Complex z)
{
   const double PI = 3.141592653589793;

   /* bf[1..N-1] are the even Bernoulli numbers / (2 n + 1)! */
   /* generated by: Table[BernoulliB[2 n]/(2 n + 1)!, {n, 1, 9}] */
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

   const double rz = creal(z);
   const double iz = cimag(z);
   const double nz = rz*rz + iz*iz;

   /* special cases */
   if (iz == 0.0) {
      if (rz <= 1.0) {
         return li2(rz);
      }
      if (rz > 1.0) {
         return li2(rz) - PI*log(rz)*I;
      }
   } else if (nz < DBL_EPSILON) {
      return z;
   }

   double _Complex cy = 0.0, cz = 0.0;
   int sgn = 1;

   /* transformation to |z|<1, Re(z)<=0.5 */
   if (rz <= 0.5) {
      if (nz > 1.0) {
         const double _Complex lz = fast_clog(-z);
         cy = -0.5*lz*lz - PI*PI/6;
         cz = -fast_clog(1.0 - 1.0 / z);
         sgn = -1;
      } else { /* nz <= 1 */
         cy = 0;
         cz = -fast_clog(1.0 - z);
         sgn = 1;
      }
   } else { /* rz > 0.5 */
      if (nz <= 2*rz) {
         cz = -fast_clog(z);
         cy = cz*fast_clog(1.0 - z) + PI*PI/6;
         sgn = -1;
      } else { /* nz > 2*rz */
         const double _Complex lz = fast_clog(-z);
         cy = -0.5*lz*lz - PI*PI/6;
         cz = -fast_clog(1.0 - 1.0 / z);
         sgn = -1;
      }
   }

   /* the dilogarithm */
   const double _Complex cz2 = cz*cz;
   const double _Complex sum =
      cz +
      cz2 * (bf[0] +
      cz  * (bf[1] +
      cz2 * (bf[2] +
      cz2 * (bf[3] +
      cz2 * (bf[4] +
      cz2 * (bf[5] +
      cz2 * (bf[6] +
      cz2 * (bf[7] +
      cz2 * (bf[8] +
      cz2 * (bf[9]))))))))));

   return (double)sgn * sum + cy;
}


/**
 * @brief Complex dilogarithm \f$\mathrm{Li}_2(z)\f$ with long double precision
 * @param z complex argument
 * @note Implementation translated from SPheno to C++
 * @return \f$\mathrm{Li}_2(z)\f$
 */
static long double _Complex cli2l(long double _Complex z)
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
   const long double nz = rz*rz + iz*iz;

   /* special cases */
   if (iz == 0.0) {
      if (rz <= 1.0) {
         return li2l(rz);
      }
      if (rz > 1.0) {
         return li2l(rz) - PI*logl(rz)*I;
      }
   } else if (nz < LDBL_EPSILON) {
      return z;
   }

   long double _Complex cy = 0.0, cz = 0.0;
   int sgn = 1;

   /* transformation to |z|<1, Re(z)<=0.5 */
   if (rz <= 0.5L) {
      if (nz > 1.0L) {
         const long double _Complex lz = fast_clogl(-z);
         cy = -0.5L*lz*lz - PI*PI/6;
         cz = -fast_clogl(1.0L - 1.0L/z);
         sgn = -1;
      } else { /* nz <= 1 */
         cy = 0;
         cz = -fast_clogl(1.0L - z);
         sgn = 1;
      }
   } else { /* rz > 0.5L */
      if (nz <= 2*rz) {
         cz = -fast_clogl(z);
         cy = cz * fast_clogl(1.0L - z) + PI*PI/6;
         sgn = -1;
      } else { /* nz > 2*rz */
         const long double _Complex lz = fast_clogl(-z);
         cy = -0.5L*lz*lz - PI*PI/6;
         cz = -fast_clogl(1.0L - 1.0L/z);
         sgn = -1;
      }
   }

   /* the dilogarithm */
   const long double _Complex cz2 = cz*cz;
   long double _Complex sum = 0.0L;

   for (int i = sizeof(bf)/sizeof(bf[0]) - 1; i >= 2; i--) {
      sum = cz2 * (bf[i] + sum);
   }

   /* lowest order terms w/ different powers */
   sum = cz + cz2 * (bf[0] + cz * (bf[1] + sum));

   return (long double)sgn * sum + cy;
}


/** C++ wrapper for cli2 */
void cli2_(double re, double im, double* res_re, double* res_im)
{
   const double _Complex result = cli2(re + I*im);
   *res_re = creal(result);
   *res_im = cimag(result);
}


/** C++ wrapper for cli2l */
void cli2l_(long double re, long double im, long double* res_re, long double* res_im)
{
   const long double _Complex result = cli2l(re + I*im);
   *res_re = creall(result);
   *res_im = cimagl(result);
}
