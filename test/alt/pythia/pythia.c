#include <math.h>

#ifndef M_PI
#define M_PI 3.1415926535897932
#endif

/*
  Dilogarithm

  Implementation from Pythia 8.315.
  Translated to C by Alexander Voigt.
 */
double pythia_dilog(const double x, const double kmax, const double xerr)
{
   if (x < 0.0) {
      return 0.5*pythia_dilog(x*x, kmax, xerr) - pythia_dilog(-x, kmax, xerr);
   }

   if (x <= 0.5) {
      double sum = x, term = x;

      for (int k = 2; k < kmax; k++) {
         double rk = (k-1.0)/k;
         term *= x*rk*rk;
         sum += term;
         if (fabs(term/sum) < xerr) return sum;
      }

      /* cout << "Maximum number of iterations exceeded in pythia_dilog" << endl; */
      return sum;
   }

   if (x < 1.0) {
      return M_PI*M_PI/6.0 - pythia_dilog(1.0 - x, kmax, xerr) - log(x)*log(1.0 - x);
   }

   if (x == 1.0) {
      return M_PI*M_PI/6.0;
   }

   if (x <= 1.01) {
      const double eps = x - 1.0,
         lne = log(eps),
         c0 = M_PI*M_PI/6.0,
         c1 =  1.0 - lne,
         c2 = -(1.0 - 2.0*lne)/4.0,
         c3 = (1.0 - 3.0*lne)/9.0,
         c4 = -(1.0 - 4.0*lne)/16.0,
         c5 = (1.0 - 5.0*lne)/25.0,
         c6 = -(1.0 - 6.0*lne)/36.0,
         c7 = (1.0 - 7.0*lne)/49.0,
         c8 = -(1.0 - 8.0*lne)/64.0;

      return c0 + eps*(c1 + eps*(c2 + eps*(c3 + eps*(c4 + eps*(c5 + eps*(c6 + eps*(c7 + eps*c8)))))));
   }

   double logx = log(x);

   if (x <= 2.0) {
      return M_PI*M_PI/6.0 + pythia_dilog(1.0 - 1.0/x, kmax, xerr) -
         logx*(log(1.0 - 1.0/x) + 0.5*logx);
   }

   return M_PI*M_PI/3.0 - pythia_dilog(1.0/x, kmax, xerr) - 0.5*logx*logx;
}
