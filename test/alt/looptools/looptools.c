#include <complex.h>
#include <math.h>


/*
  Implementation of dilogarithm from FeynHiggs 2.16.0.

  file: src/LT/spence.F

  Translated to C by Alexander Voigt
 */
static long double complex spence(long double complex x)
{
   return x;
}


void looptools_dilog(long double re, long double im, long double* res_re, long double* res_im)
{
   long double complex z = re + I*im;
   long double complex result = spence(z);
   *res_re = creall(result);
   *res_im = cimagl(result);
}
