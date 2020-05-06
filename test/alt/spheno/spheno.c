#include <complex.h>
#include <float.h>
#include <math.h>

/* --------------------------------------------------------------------- */
/*  Spence-fkt. li2:                                                     */
/*  local double prec. Version of function dilogarithm, 22. 5. 95        */
/*  H. Eberl, Convention - see Mathematica-manual side 575               */
/*  checked                                                              */
/*  1.10.2000: change to f90 + implicit none statement                   */
/* --------------------------------------------------------------------- */
static double Li2(double x)
{
   const double PI  = 3.141592653589793;
   const double HF  = 0.5;
   const double PI2 = PI*PI;
   const double PI3 = PI2/3;
   const double PI6 = PI2/6;
   const double PI12 = PI2/12;
   const double C[20] = {
       0.42996693560813697,
       0.40975987533077105,
      -0.01858843665014592,
       0.00145751084062268,
      -0.00014304184442340,
       0.00001588415541880,
      -0.00000190784959387,
       0.00000024195180854,
      -0.00000003193341274,
       0.00000000434545063,
      -0.00000000060578480,
       0.00000000008612098,
      -0.00000000001244332,
       0.00000000000182256,
      -0.00000000000027007,
       0.00000000000004042,
      -0.00000000000000610,
       0.00000000000000093,
      -0.00000000000000014,
       0.00000000000000002
   };

   double H,T,Y,S,A,ALFA,B0,B1,B2;

   if (x == 1) {
       H = PI6;
   } else if (x == -1) {
       H = -PI12;
   } else {
       T = -x;
       if (T <= -2) {
           Y = -1/(1+T);
           S = 1;
           B1= log(-T);
           B2= log(1+1/T);
           A = -PI3+HF*(B1*B1-B2*B2);
       } else if (T < -1) {
           Y = -1-T;
           S = -1;
           A = log(-T);
           A = -PI6+A*(A+log(1+1/T));
       } else if (T <= -0.5) {
           Y = -(1+T)/T;
           S = 1;
           A = log(-T);
           A = -PI6+A*(-HF*A+log(1+T));
       } else if (T < 0) {
           Y = -T/(1+T);
           S = -1;
           B1= log(1+T);
           A = HF*B1*B1;
       } else if (T <= 1) {
           Y = T;
           S = 1;
           A = 0;
       } else {
           Y = 1/T;
           S = -1;
           B1= log(T);
           A = PI6+HF*B1*B1;
       }
       H    = Y+Y-1;
       ALFA = H+H;
       B1   = 0;
       B2   = 0;
       for (int i = 19; i >= 0; i--) {
          B0 = C[i] + ALFA*B1-B2;
          B2 = B1;
          B1 = B0;
       }
       H = -(S*(B0-H*B2)+A);
    }
    return H;
}

/*
  Implementation of dilogarithm from SPheno 4.0.4
  Written by Werner Porod

  Translated to C by Alexander Voigt
 */
static double complex CLI2(double complex z)
{
   const double Pi = 3.141592653589793;
   const double pi2o12 = Pi * Pi / 12.0;

   /* --------------------------------------------- */
   /*  Bernulli numbers / (2 n + 1)!                */
   /* --------------------------------------------- */
   const double bf[20] = {
      -1. / 4.,
      +1. / 36.,
      -1. / 36.e2,
      +1. / 21168.e1,
      -1. / 108864.e2,
      +1. / 52690176.e1,
      -691. / 16999766784.e3,
      +1. / 1120863744.e3,
      -3617. / 18140058832896.e4,
      +43867. / 97072790126247936.e3,
      -174611. / 168600109166641152.e5,
      +77683. / 32432530090601152512.e4,
      -236364091. / 4234560341829359173632.e7,
      +657931. / 5025632054039239458816.e6,
      -3392780147. / 109890470493622010006470656.e7,
      +172.3168255201 / 2355349904102724211909.3102313472e6,
      -770.9321041217 / 4428491985594062112714.2791446528e8,
      4.1576356446138997196178996207752266734882541595116e-29,
      -9.9621484882846221031940067024558388498548600173945e-31,
      2.3940344248961653005211679878937495629342791569329e-32
   };

   const double rz = creal(z);
   const double iz = cimag(z);
   const double az = sqrt(rz*rz + iz*iz);

   /* ------------------ */
   /*  special cases     */
   /* ------------------ */
   if (iz == 0.) {
      if (rz <= 1.0)
         return Li2(rz) + I*0.;
      if (rz > 1.0)
         return Li2(rz) - I*Pi*log(rz);
   } else if (az < DBL_EPSILON) {
      return z;
   }

   double complex cy, cz;
   int jsgn, ipi12;

   /* ----------------------------------------------------- */
   /*  transformation to |z|<1, Re(z)<=0.5                  */
   /* ----------------------------------------------------- */
   if (rz <= 0.5) {
      if (az > 1.0) {
         double complex l = clog(-z);
         cy = -0.5 * l*l;
         cz = -clog(1.0 - 1.0 / z);
         jsgn = -1;
         ipi12 = -2;
      } else { /* az <= 1.0 */
         cy = 0;
         cz = -clog(1.0 - z);
         jsgn = 1;
         ipi12 = 0;
      }
   } else { /* rz > 0.5 */
      if (az <= sqrt(2*rz)) {
         cz = -clog(z);
         cy = cz * clog(1.0 - z);
         jsgn = -1;
         ipi12 = 2;
      } else { /* az > sqrt(2*rz) */
         double complex l = clog(-z);
         cy = -0.5 * l*l;
         cz = -clog(1.0 - 1.0 / z);
         jsgn = -1;
         ipi12 = -2;
      }
   }

   /* -------------------------------------------- */
   /*  the dilogarithm                             */
   /* -------------------------------------------- */
   const double complex cz2 = cz*cz;
   double complex sumC = cz2 * bf[19];

   for (int i1 = 4; i1 <= 20; i1++)
      sumC = cz2 * (sumC + bf[22 - i1]);

   /* watch the powers of z */
   sumC = cz + cz2 * (bf[0] + cz * (bf[1] + sumC));

   complex double CLi2;

   if (jsgn == 1)
      CLi2 = sumC + cy;
   else
      CLi2 = -sumC + cy;

   if (ipi12 != 0)
      CLi2 = CLi2 + ipi12 * pi2o12;

   return CLi2;
}


void spheno_dilog(double re, double im, double* res_re, double* res_im)
{
   double complex z = re + I*im;
   double complex result = CLI2(z);
   *res_re = creal(result);
   *res_im = cimag(result);
}
