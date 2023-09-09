/* ====================================================================
 * This file is part of Polylogarithm.
 *
 * Polylogarithm is licenced under the MIT License.
 * ==================================================================== */

#include "c_wrappers.h"
#include <complex.h>

float _Complex cli2f(float _Complex);

double _Complex cli2(double _Complex);
double _Complex cli3(double _Complex);
double _Complex cli4(double _Complex);
double _Complex cli5(double _Complex);
double _Complex cli6(double _Complex);

long double _Complex cli2l(long double _Complex);
long double _Complex cli3l(long double _Complex);
long double _Complex cli4l(long double _Complex);
long double _Complex cli5l(long double _Complex);
long double _Complex cli6l(long double _Complex);

/** C++ wrapper for cli2f */
void cli2f_c(float re, float im, float* res_re, float* res_im)
{
   const float _Complex result = cli2f(re + I*im);
   *res_re = crealf(result);
   *res_im = cimagf(result);
}

/** C++ wrapper for cli2 */
void cli2_c(double re, double im, double* res_re, double* res_im)
{
   const double _Complex result = cli2(re + I*im);
   *res_re = creal(result);
   *res_im = cimag(result);
}

/** C++ wrapper for cli2l */
void cli2l_c(long double re, long double im, long double* res_re, long double* res_im)
{
   const long double _Complex result = cli2l(re + I*im);
   *res_re = creall(result);
   *res_im = cimagl(result);
}

/** C++ wrapper for cli3 */
void cli3_c(double re, double im, double* res_re, double* res_im)
{
   const double _Complex result = cli3(re + I*im);
   *res_re = creal(result);
   *res_im = cimag(result);
}

/** C++ wrapper for cli3l */
void cli3l_c(long double re, long double im, long double* res_re, long double* res_im)
{
   const long double _Complex result = cli3l(re + I*im);
   *res_re = creall(result);
   *res_im = cimagl(result);
}

/** C++ wrapper for cli4 */
void cli4_c(double re, double im, double* res_re, double* res_im)
{
   const double _Complex result = cli4(re + I*im);
   *res_re = creal(result);
   *res_im = cimag(result);
}

/** C++ wrapper for cli4l */
void cli4l_c(long double re, long double im, long double* res_re, long double* res_im)
{
   const long double _Complex result = cli4l(re + I*im);
   *res_re = creall(result);
   *res_im = cimagl(result);
}

/** C++ wrapper for cli5 */
void cli5_c(double re, double im, double* res_re, double* res_im)
{
   const double _Complex result = cli5(re + I*im);
   *res_re = creal(result);
   *res_im = cimag(result);
}

/** C++ wrapper for cli5l */
void cli5l_c(long double re, long double im, long double* res_re, long double* res_im)
{
   const long double _Complex result = cli5l(re + I*im);
   *res_re = creall(result);
   *res_im = cimagl(result);
}

/** C++ wrapper for cli6 */
void cli6_c(double re, double im, double* res_re, double* res_im)
{
   const double _Complex result = cli6(re + I*im);
   *res_re = creal(result);
   *res_im = cimag(result);
}

/** C++ wrapper for cli6l */
void cli6l_c(long double re, long double im, long double* res_re, long double* res_im)
{
   const long double _Complex result = cli6l(re + I*im);
   *res_re = creall(result);
   *res_im = cimagl(result);
}
