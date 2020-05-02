#include <complex.h>
#include <math.h>


/* transforms -0.0 -> 0.0 */
static double complex pclog(double complex z)
{
   double re = creal(z) == 0.0 ? 0.0 : creal(z);
   double im = cimag(z) == 0.0 ? 0.0 : cimag(z);
   return clog(re + I*im);
}

/*
  Implementation of dilogarithm from Sherpa 2.2.9

  file: Examples/H_in_TTBar/LHC_TTH_MCatNLO/lib.f
  file: Examples/V_plus_Bs/LHC_Wbb/Wbb_Virtual.f

  Translated to C by Alexander Voigt
 */
static double complex spence(double complex x)
{
   double pi, xr, xi, vor;
   complex double z, d, p, result;
   const double a[19] = {
      [0]  = 1.0,
      [1]  = -0.5,
      [2]  = 1.0/6.0,
      [4]  = -1.0/30.0,
      [6]  = 1.0/42.0,
      [8]  = -1.0/30.0,
      [10] = 5.0/66.0,
      [12] = -691.0/2730.0,
      [14] = 7.0/6.0,
      [16] = -3617.0/510.0,
      [18] = 43867.0/798.0,
   };

   pi = 4.0*atan(1.0);
   xr = creal(x);
   xi = cimag(x);

   if (xr != 1.0) goto L111;
   if (xi == 0.0) goto L20;

L111:
   /* projection into the convergence radius */
   vor = 1.0;
   p = 0.0 + I*0.0;

   if (creal(x) <= 0.5) goto L1;

   p = pi*pi/6.0 - pclog(x)*pclog(1.0 - x);
   vor = -1.0;
   x = 1.0 - x;

L1:
   if (cabs(x) < 1.0) goto L2;

   p = p - (pi*pi/6.0 + pclog(-x)*pclog(-x)/2.0)*vor;
   vor = vor*(-1.0);
   x = 1.0/x;

L2:
   /* calculation of the spence function */
   z = (-1.0)*pclog(1.0 - x);
   d = a[18] + I*0.0;

   for (int n = 1; n <= 18; n++) {
      d = d*z/(20.0 - n) + a[18 - n];
   }

   d = d*z;
   result = d*vor + p;
   goto L30;

L20:
   result = pi*pi/6.0;

L30:
   return result;
}


void sherpa_dilog(double re, double im, double* res_re, double* res_im)
{
   double complex z = re + I*im;
   double complex result = spence(z);
   *res_re = creal(result);
   *res_im = cimag(result);
}
