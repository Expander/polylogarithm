#include <tsil.h>

#ifdef ENABLE_FORTRAN

void dli2_wrapper(const double*, double*);
void cdli2_wrapper(const double*, const double*, double*, double*);
void cdli3_wrapper(const double*, const double*, double*, double*);

double li2_fortran(double x)
{
   double res = 0;
   dli2_wrapper(&x, &res);
   return res;
}

void cli2_fortran(double re, double im, double* res_re, double* res_im)
{
   cdli2_wrapper(&re, &im, res_re, res_im);
}

void cli3_fortran(double re, double im, double* res_re, double* res_im)
{
   cdli3_wrapper(&re, &im, res_re, res_im);
}

#endif

TSIL_REAL TSIL_dilog_real(TSIL_REAL x);

long double tsil_dilog_real(long double x)
{
   return TSIL_dilog_real(x);
}

void tsil_dilog_complex(long double re, long double im, long double* res_re, long double* res_im)
{
   TSIL_COMPLEX z = re + I*im;
   TSIL_COMPLEX res = TSIL_Dilog(z);
   *res_re = creall(res);
   *res_im = cimagl(res);
}

void tsil_trilog_complex(long double re, long double im, long double* res_re, long double* res_im)
{
   TSIL_COMPLEX z = re + I*im;
   TSIL_COMPLEX res = TSIL_Trilog(z);
   *res_re = creall(res);
   *res_im = cimagl(res);
}
