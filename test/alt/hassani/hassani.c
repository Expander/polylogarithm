#include <math.h>

static const double PI = 3.141592653589793;

/*
  Series expansion of Spence function for real x > 1

  Mehdi Hassani - Approximation of the dilogarithm function.
  INSTITUTE FOR ADVANCED STUDIES IN BASIC SCIENCES
  P.O. B OX 45195-1159
  ZANJAN, IRAN.
  Received 15 April, 2006; accepted 03 January, 2007
  Communicated by A. Lupaş

  Translated to C by Alexander Voigt
 */
static double hassani_spence(double x)
{
   const double lx = log(x);

   double sum = 0.0;
   double xn = 1;
   const int N = 50; // Log[1/2, 10^-15]

   for (int i = 1; i < N; i++) {
      xn *= x;
      sum += (1.0/i + lx)/(i*xn);
   }

   return -0.5*lx*lx - PI*PI/6 + sum;
}

double hassani_dilog(double x)
{
   /* tranform to x < 0 */
   if (x > 2) {
      return -hassani_dilog(1 - x) + PI*PI/6 - log(x)*log(x - 1);
   } else if (x > 1) {
      const double l = log(x - 1);
      return hassani_dilog(1/(1 - x)) + PI*PI/3 + l*(0.5*l - log(x));
   } else if (x == 1) {
      return PI*PI/6;
   } else if (x > 0.5) {
      const double l = log(1 - x);
      return -hassani_dilog(x/(x - 1)) - 0.5*l*l;
   } else if (x == 0.5) {
      const double l2 = 0.69314718055994531;
      return PI*PI/12 - 0.5*l2*l2;
   } else if (x > 0) {
      const double l = log(x);
      return hassani_dilog(1 - 1/x) + PI*PI/6 - l*(log(1/x - 1) + 0.5*l);
   } else if (x == 0) {
      return 0;
   } else if (x > -1) {
      const double l = log(-x);
      return -hassani_dilog(1/x) - PI*PI/6 - 0.5*l*l;
   }

   return hassani_spence(1 - x);
}
