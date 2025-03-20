#include <math.h>

/*
  Algorithm 490
  The Dilogarithm Function of a Real Argument [S22]
  Edward S. Ginsberg (Recd 22 June 1973)
  Department of Physics, University of Massachusetts
  at Boston, Boston, MA 02125
  and
  Dorothy Zaborowski
  Information Processing Center, Massachusetts Institute
  of Technology, Cambridge, MA 02139

  Published in Communications of the ACM,
  Volume 18, Number 4, p. 200f, April 1975

  Translated to C by Alexander Voigt
 */
double algorithm_490(double X)
{
   const double PI  = 3.141592653589793;
   const double PI2 = PI*PI;
   const double PI3 = PI2/3;
   const double PI6 = PI2/6;
   const double C[30] = {
      36,        576,       3600,      14400,     44100,     112896,
      254016,    518400,    980100,    1742400,   2944656,   4769856,
      7452900,   11289600,  16646400,  23970816,  33802596,  46785600,
      63680400,  85377600,  112911876, 147476736, 190440000, 243360000,
      308002500, 386358336, 480661776, 593409600, 727380900, 885657600
   };

   double A, B, BY, C1, C2, C3, C4, DX, DY, TEST, W, X0, Y, Z, DILOG;

   if (X > 12.6) goto L10;
   if (X >= 12.59) goto L100;
   if (X >= 2) goto L10;
   if (X > 1) goto L20;
   if (X == 1) goto L30;
   if (X > 0.5) goto L40;
   if (X > 0.01) goto L50;
   if (X < -1) goto L60;
   if (X < -0.01) goto L70;

   DILOG = X*(1 + X*(0.25 + X*(1.0/9 + X*(0.0625 + X*(0.04 + X*(1.0/36 + X*(1.0/49 + X/64)))))));
   return DILOG;

L10:
   Y = 1/X;
   BY = -1 - Y*(4 + Y);
   DX = log(X);
   DILOG = PI3 - 0.5*DX*DX + (Y*(4 + 5.75*Y) + 3*(1 + Y)*(1 - Y)*log(1 - Y))/BY;
   if (DILOG + 4*Y == DILOG)
      return DILOG;
   goto L80;

L20:
   Y = 1 - 1/X;
   DX = log(X);
   BY = 1 + Y*(4 + Y);
   DILOG = PI6 + DX*(0.5*DX - log(X - 1)) + (Y*(4 + 5.75*Y) - 3*(1 + Y)*DX/X)/BY;
   goto L80;

L30:
   DILOG = PI6;
   return DILOG;

L40:
   Y = 1 - X;
   DX = log(X);
   BY = -1 - Y*(4 + Y);
   DILOG = PI6 - DX*log(Y) + (Y*(4 + 5.75*Y) + 3*(1 + Y)*DX*X)/BY;
   goto L80;

L50:
   Y = X;
   BY = 1 + Y*(4 + Y);
   DILOG = (Y*(4 + 5.75*Y) + 3*(1 + Y)*(1 - Y)*log(1 - Y))/BY;
   goto L80;

L60:
   Y = 1/(1 - X);
   DX = log(-X);
   DY = log(Y);
   BY = 1 + Y*(4 + Y);
   DILOG = -PI6 + 0.5*DY*(DY + 2*DX) + (Y*(4 + 5.75*Y) + 3*(1 + Y)*(1 - Y)*(DX + DY))/BY;
   if (DILOG + 4*Y == DILOG)
      return DILOG;
   goto L80;

L70:
   Y = X/(X - 1);
   DX = log(1 - X);
   BY = -1 - Y*(4 + Y);
   DILOG = (Y*(4 + 5.75*Y) - 3*(1 + Y)*(1 - Y)*DX)/BY - 0.5*DX*DX;

L80:
   B = 4*Y*Y/BY;
   for (int N = 0; N < 30; N++) {
      B = B*Y;
      A = B/C[N];
      TEST = DILOG;
      DILOG = DILOG + A;
      if (DILOG == TEST)
         return DILOG;
   }
   return DILOG;

L100:
   X0 = 12.5951703698450161286398965;
   Y = X/X0 - 1;
   Z = 1/11.5951703698450161286398965;
   W = Y*Z;
   C1 = (3*X0 - 2)/6;
   C2 = ((11*X0 - 15)*X0 + 6)/24;
   C3 = (((50*X0 - 104)*X0 + 84)*X0 - 24)/120;
   C4 = ((((274*X0 - 770)*X0 + 940)*X0 - 540)*X0 + 120)/720;
   DILOG = Y*(1 - Y*(0.5 - Y*(1.0/3.0 - Y*(0.25 - Y*(0.2 - Y/6)))))*log(Z) - W*X0*Y*(0.5 - W*(C1 - W*(C2 - W*(C3 - W*C4))));

   return DILOG;
}


/*
 * Refined version of algorithm 490 without goto statements and with
 * one division less.
 *
 * @author Refined by Alexander Voigt
 */
double algorithm_490_2(double x)
{
   const double PI = 3.1415926535897932;
   // Table[1/(n(n+1)(n+2))^2, {n,2,25}] // N[#,17]& // CForm
   const double C[24] = {
      0.0017361111111111111, 0.00027777777777777778,
      0.000069444444444444444, 0.000022675736961451247,
      8.8577097505668934e-6, 3.9367598891408415e-6,
      1.9290123456790123e-6, 1.0203040506070809e-6,
      5.7392102846648301e-7, 3.3959824169614379e-7,
      2.0964993492466020e-7, 1.3417595835178253e-7,
      8.8577097505668934e-8, 6.0073048827374087e-8,
      4.1717395019009783e-8, 2.9583526661680067e-8,
      2.1374098013063849e-8, 1.5703418948373440e-8,
      1.1712674050336388e-8, 8.8564643102732612e-9,
      6.7807304875529656e-9, 5.2509976895610166e-9,
      4.1091387245233399e-9, 3.2467268934505402e-9
   };

   double y = 0, r = 0, s = 1;

   /* transform to [0, 1/2] */
   if (x < -1) {
      const double l = log(1 - x);
      y = 1/(1 - x);
      r = -PI*PI/6 + l*(0.5*l - log(-x));
      s = 1;
   } else if (x == -1) {
      return -PI*PI/12;
   } else if (x < 0) {
      const double l = log1p(-x);
      y = x/(x - 1);
      r = -0.5*l*l;
      s = -1;
   } else if (x == 0) {
      return x;
   } else if (x < 0.5) {
      y = x;
      r = 0;
      s = 1;
   } else if (x < 1) {
      y = 1 - x;
      r = PI*PI/6 - log(x)*log1p(-x);
      s = -1;
   } else if (x == 1) {
      return PI*PI/6;
   } else if (x < 2) {
      const double l = log(x);
      y = 1 - 1/x;
      r = PI*PI/6 - l*(log(y) + 0.5*l);
      s = 1;
   } else {
      const double l = log(x);
      y = 1/x;
      r = PI*PI/3 - 0.5*l*l;
      s = -1;
   }

   double yn = 4*y*y*y;
   double sum = yn/36;

   for (int n = 0; n < 24; ++n) {
      const int d = (n + 2)*(n + 3)*(n + 4);
      yn *= y;
      sum += yn/(d*d);
   }

   const double d = 1 + y*(4 + y);

   return r + s*(sum + y*(4 + 5.75*y) + 3*(1 + y)*(1 - y)*log(1 - y))/d;
}
