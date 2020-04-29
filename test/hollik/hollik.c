#include <math.h>
#include <complex.h>

/*
  SPENCE-FUNKTION KOMPLEX, FREI NACH HOLLIK
  20.07.83
  LAST CHANGED 10.05.89, ANSGAR DENNER

  Translated to C by Alexander Voigt
 */
double complex pCSPEN(double complex Z)
{
   double complex W,SUM,U,L;
   double RZ,AZ,A1;

   // BEACHTE:                 B(N)=B2N
   // B(1)=1./6.
   // B(2)=-1./30.
   // B(3)=1./42.
   // B(4)=-1./30.
   // B(5)=5./66.
   // B(6)=-691./2730.
   // B(7)=7./6.
   // B(8)=-3617./510.
   // B(9)=43867./798.
   // PI=3.1415926535897932384
   // PI*PI/6.=1.6449..., PI*PI/3=3.28986...
   const double B[9] = {
       0.1666666666666666666666666667,
      -0.0333333333333333333333333333,
       0.0238095238095238095238095238,
      -0.0333333333333333333333333333,
       0.0757575757575757575757575758,
      -0.2531135531135531135531135531,
       1.1666666666666666666666666667,
      -7.09215686274509804,
      54.97117794486215539
   };

   RZ = creal(Z);
   AZ = cabs(Z);
   A1 = cabs(1-Z);

   if (AZ < 1E-20)
      return -clog(1-Z);

   if ((RZ == 1.0) && (cimag(Z) == 0.0))
      return 1.64493406684822643;

   if (RZ > 5E-1) goto L20;
   if (AZ > 1) goto L10;

   W = -clog(1-Z);
   SUM = W - 0.25*W*W;
   U = W;

   if (cabs(U) < 1E-10) goto L2;

   for (int K = 1; K <= 9; K++) {
      U = U*W*W/(2*K*(2*K+1));
      if (cabs(U*B[K-1]/SUM) < 1E-20)
         break;
      SUM = SUM + U*B[K-1];
   }

 L2:
   return SUM;

 L10:
   W = -clog(1-1/Z);
   SUM = W - 0.25*W*W;
   U = W;

   if (cabs(U) < 1E-10) goto L12;

   for (int K = 1; K <= 9; K++) {
      U = U*W*W/(2*K*(2*K+1));
      if (cabs(B[K-1]*U/SUM) < 1E-20)
         break;
      SUM = SUM + U*B[K-1];
   }
 L12:
   L = clog(-Z);

   return -SUM - 1.64493406684822643 - 0.5*L*L;

 L20:
   if (A1 > 1) goto L30;

   W = -clog(Z);
   SUM = W - 0.25*W*W;
   U = W;

   if (cabs(U) < 1E-10) goto L22;

   for (int K = 1; K <= 9; K++) {
      U = U*W*W/(2*K*(2*K+1));
      if (cabs(U*B[K-1]/SUM) < 1E-20)
         break;
      SUM = SUM + U*B[K-1];
   }
 L22:
   return -SUM + 1.64493406684822643 - clog(Z)*clog(1-Z);

 L30:
   W = clog(1-1/Z);
   SUM = W - 0.25*W*W;
   U = W;

   if (cabs(U) < 1E-10) goto L32;

   for (int K = 1; K <= 9; K++) {
      U = U*W*W/(2*K*(2*K+1));
      if (cabs(U*B[K-1]/SUM) < 1E-20)
         break;
      SUM = SUM + U*B[K-1];
   }
 L32:
   L = clog(Z-1);

   return SUM + 3.28986813369645287
      + 0.5*L*L - clog(Z)*clog(1-Z);
}


void hollik_dilog(double re, double im, double* res_re, double* res_im)
{
   double complex z = re + I*im;
   double complex result = pCSPEN(z);
   *res_re = creal(result);
   *res_im = cimag(result);
}
