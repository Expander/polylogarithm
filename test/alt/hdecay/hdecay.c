#include <math.h>
#include <complex.h>


/*
  Complex dilogarithm from HDECAY
  Authors: A. Djouadi, J. Kalinowski, M. Spira
  License: GPL 2.0

  Translated to C by Alexander Voigt
 */
double complex _hdecay_dilog(double complex z)
{
   return z;
}



void hdecay_dilog(double re, double im, double* res_re, double* res_im)
{
   double complex z = re + I*im;
   double complex result = _hdecay_dilog(z);
   *res_re = creal(result);
   *res_im = cimag(result);
}
