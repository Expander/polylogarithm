#include <math.h>

/**
 * @brief Clausen function \f$\operatorname{Cl}_2(\theta)\f$
 * @param x real angle
 * @return \f$\operatorname{Cl}_2(\theta)\f$
 * @author Alexander Voigt
 *
 * Implementation as series expansion from Jiming Wu, Xiaoping Zhang,
 * Dongjie Liu: "An efficient calculation of the Clausen functions
 * Cl_n(0)(n >= 2)"
 */
double clausen_2_wu(double x)
{
   const double PI = 3.14159265358979324;
   const double PI2 = 2*PI;
   const double C[] = {
      1.0,
      -2.7777777777777778e-02,
      -2.7777777777777778e-04,
      -4.7241118669690098e-06,
      -9.1857730746619636e-08,
      -1.8978869988970999e-09,
      -4.0647616451442255e-11,
      -8.9216910204564526e-13,
      -1.9939295860721076e-14,
      -4.5189800296199182e-16,
      -1.0356517612181247e-17,
      -2.3952186210261867e-19,
      -5.5817858743250093e-21,
      -1.3091507554183213e-22,
      -3.0874198024267403e-24,
      -7.3159756527022034e-26,
      -1.7408456572340007e-27,
      -4.1576356446138997e-29,
      -9.9621484882846221e-31,
      -2.3940344248961653e-32,
      -5.7683473553673901e-34,
      -1.3931794796470080e-35,
      -3.3721219654850895e-37,
      -8.1782087775621026e-39,
      -1.9870108311523859e-40
   };

   double sgn = 1;

   if (x < 0) {
      x = -x;
      sgn = -1;
   }

   if (x >= PI2) {
      x = fmod(x, PI2);
   }

   if (x > PI) {
      const double p0 = 6.28125;
      const double p1 = 0.0019353071795864769253;
      x = (p0 - x) + p1;
      sgn = -sgn;
   }

   if (x == 0 || x == PI) {
      return 0;
   }

   const double x2 = x*x;
   double sum = C[sizeof(C)/sizeof(C[0]) - 1];

   for (int i = sizeof(C)/sizeof(C[0]) - 2; i >= 0; --i) {
      sum = x2*sum + C[i];
   }

   sum -= log(2*sin(x/2));

   return sgn*x*sum;
}

/**
 * @brief Clausen function \f$\operatorname{Cl}_3(\theta)\f$
 * @param x real angle
 * @return \f$\operatorname{Cl}_3(\theta)\f$
 * @author Alexander Voigt
 *
 * Implementation as series expansion from Jiming Wu, Xiaoping Zhang,
 * Dongjie Liu: "An efficient calculation of the Clausen functions
 * Cl_n(0)(n >= 2)"
 */
double clausen_3_wu(double x)
{
   const double PI = 3.14159265358979324;
   const double PI2 = 2*PI;
   const double C[] = {
      1.2020569031595943,
     -0.75,
      1.7361111111111111e-02,
      1.6203703703703704e-04,
      2.6573129251700680e-06,
      5.0521751910640800e-08,
      1.0280221244025958e-09,
      2.1775508813272637e-11,
      4.7396483546174904e-13,
      1.0523517259825012e-14,
      2.3724645155504571e-16,
      5.4136342063674700e-18,
      1.2475096984511389e-19,
      2.8982349732072164e-21,
      6.7795306977020209e-23,
      1.5951668979204825e-24,
      3.7722999459245736e-26,
      8.9602350004691215e-28,
      2.1365627618154762e-29,
      5.1121551453039508e-31,
      1.2269426427592847e-32,
      2.9528444795333068e-34,
      7.1242132481949272e-36,
      1.7227144823673827e-37,
      4.1742940635473232e-39
   };

   if (x < 0) {
      x = -x;
   }

   if (x >= PI2) {
      x = fmod(x, PI2);
   }

   if (x > PI) {
      const double p0 = 6.28125;
      const double p1 = 0.0019353071795864769253;
      x = (p0 - x) + p1;
   }

   if (x == 0) {
      return C[0];
   }

   if (x == PI) {
      return -0.75*C[0];
   }

   const double x2 = x*x;
   double sum = C[sizeof(C)/sizeof(C[0]) - 1];

   for (int i = sizeof(C)/sizeof(C[0]) - 2; i >= 0; --i) {
      sum = x2*sum + C[i];
   }

   sum += x2/2*log(2*sin(x/2));

   return sum;
}
