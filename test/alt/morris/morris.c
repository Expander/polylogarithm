#include <math.h>


static double horner(double x, const double* c, int len)
{
   double p = 0.0;
   while (len--)
      p = p*x + c[len];
   return p;
}


/*
  Real dilogarithm function

  Based on Robert Morris,
  "The Dilogarithm Function of a Real Argument"
  MATHEMATICS OF COMPUTATION, VOLUME 33, NUMBER 146
  APRIL 1979, PAGES 778-787,
  DILOG 0011

  Translated to C by Alexander Voigt
*/
double morris_dilog(double x)
{
   const double PI = 3.141592653589793;
   const double cp[] = {
      +0.12050991261341671705797e4,
      -0.32407609254539474989552e4,
      +0.31928976938241289723624e4,
      -0.13920671465266965241768e4,
      +0.25212662919363406616773e3,
      -0.13120034432716341584577e2,
   };
   const double cq[] = {
      +0.12050991261341672354682e4,
      -0.35420357069875166000880e4,
      +0.39445067176691511608434e4,
      -0.20599529983831116588803e4,
      +0.50200962202768116987420e3,
      -0.48063997258736084391455e2,
      +0.10000000000000000000000e1,
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
      const double l = log(1 - x);
      y = 1/(1 - x);
      r = -PI*PI/6 + l*(0.5*l - log(-x));
      s = 1;
   } else if (x < 0) {
      const double l = log1p(-x);
      y = x/(x - 1);
      r = -0.5*l*l;
      s = -1;
   } else if (x < 0.5) {
      y = x;
      r = 0;
      s = 1;
   } else if (x < 1) {
      y = 1 - x;
      r = PI*PI/6 - log(x)*log(1 - x);
      s = -1;
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

   const double p = horner(y, cp, sizeof(cp)/sizeof(cp[0]));
   const double q = horner(y, cq, sizeof(cq)/sizeof(cq[0]));

   return r + s*y*p/q;
}
