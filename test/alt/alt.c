#include <tsil.h>

void li2sub_(double*, const double*);
void ltini_(void);

void init_looptools()
{
   ltini_();
}

double looptools_dilog(double x)
{
   double res = 0;
   li2sub_(&res, &x);
   return res;
}

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
