#include <complex.h>
#include <math.h>


/*
  Implementation of dilogarithm from FeynHiggs 2.16.0.

  file: src/LT/spence.F

  Translated to C by Alexander Voigt
 */
static double complex spence(double complex x)
{
   return x;
}


void looptools_dilog(double re, double im, double* res_re, double* res_im)
{
   double complex z = re + I*im;
   double complex result = spence(z);
   *res_re = creal(result);
   *res_im = cimag(result);
}
