#include <math.h>
#include <complex.h>


/* transforms -0.0 -> 0.0 */
static double complex clog_hdec(double complex z)
{
   double re = creal(z) == 0.0 ? fabs(creal(z)) : creal(z);
   double im = cimag(z) == 0.0 ? fabs(cimag(z)) : cimag(z);
   return clog(re + I*im);
}


static double complex sqr_hdec(double complex z)
{
   return z*z;
}


/* DOUBLE PRECISION VERSION OF FACTORIAL */
double factrl_hdec(int n)
{
   double result = 1.0;

   if (n == 0) {
      return result;
   }

   for (int i = 1; i <= n; i++) {
      result = result * (double)i;
   }

   return result;
}


/* TAYLOR-EXPANSION FOR COMPLEX DILOGARITHM (SPENCE-FUNCTION) */
static double complex cli2_hdec(double complex x)
{
   const double b[18] = {
     -1.0/2.0,
      1.0/6.0,
      0.0,
     -1.0/30.0,
      0.0,
      1.0/42.0,
      0.0,
     -1.0/30.0,
      0.0,
      5.0/66.0,
      0.0,
     -691.0/2730.0,
      0.0,
      7.0/6.0,
      0.0,
     -3617.0/510.0,
      0.0,
      43867.0/798.0
   };

   double b2[18];

   for (int i = 0; i < 18; i++) {
      b2[i] = b[i]/factrl_hdec(i + 2);
   }

   const int nber = 18;
   const int n = nber - 1;
   const double complex z = -clog_hdec(1.0 - x);
   double complex result = b2[nber - 1];

   for (int i = n - 1; i >= 0; i--) {
      result = z*result + b2[i];
   }

   result = sqr_hdec(z)*result + z;

   return result;
}


/*
  Complex dilogarithm from HDECAY
  Authors: A. Djouadi, J. Kalinowski, M. Spira
  License: GPL 2.0

  Translated to C by Alexander Voigt
 */
double complex li2_hdec(double complex x)
{
   const double pi = 4.0*atan(1.0);
   const double zeta2 = pi*pi/6;
   const double zero = 1e-16;
   const double xr = creal(x);
   const double xi = cimag(x);
   const double r2 = xr*xr + xi*xi;

   if (r2 < zero) {
      return x;
   }

   const double rr = xr/r2;
   double complex y = 0.0;

   if (r2 == 1.0 && xi == 0.0) {
      if (xr == 1.0) {
         return zeta2;
      } else {
         return -zeta2/2.0;
      }
   } else if (r2 > 1.0 && rr > 0.5) {
      y = (x - 1.0)/x;
      return cli2_hdec(y) + zeta2 - clog_hdec(x)*clog_hdec(1.0 - x) + 0.5*sqr_hdec(clog_hdec(x));
   } else if (r2 > 1.0 && rr <= 0.5) {
      y = 1.0/x;
      return -cli2_hdec(y) - zeta2 - 0.5*sqr_hdec(clog_hdec(-x));
   } else if (r2 <= 1.0 && xr > 0.5) {
      y = 1.0 - x;
      return -cli2_hdec(y) + zeta2 - clog_hdec(x)*clog_hdec(1.0 - x);
   } else if (r2 <= 1.0 && xr <= 0.5) {
      y = x;
      return cli2_hdec(y);
   }

   return 0.0;
}


void hdecay_dilog(double re, double im, double* res_re, double* res_im)
{
   double complex z = re + I*im;
   double complex result = li2_hdec(z);
   *res_re = creal(result);
   *res_im = cimag(result);
}
