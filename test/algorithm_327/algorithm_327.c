#include <math.h>

/*
  Algorithm 327, Dilogarithm [S22]
  K.S. Kölbig (Rect. 10 Oct. 1967)
  Applied Mathematics Group, Data Handling Division,
  European Organization for Nuclear Research (CERN)
  1211 Genea 23, Switzerland

  Published in Communications of the ACM,
  Volume 11, Number 4, p. 270f, April 1968

  Translated to C by Alexander Voigt
 */
double algorithm_327(double x)
{
   const double PI  = 3.141592653589793;
   const double PI2 = PI*PI;
   const double PI3 = PI2/3;
   const double PI6 = PI2/6;

   double f, u, y, z, l;

   if (x >= 2) {
      l = log(x);
      z = 1/x;
      u = -0.5*l*l + PI3;
      f = -1;
   } else if (x > 1) {
      z = (x - 1)/x;
      u = -0.5*log(x)*log(z*x - z) + PI6;
      f = 1;
   } else if (x == 1) {
      return PI6;
   } else if (x > 0.5) {
      z = 1 - x;
      u = -log(x)*log(z) + PI6;
      f = -1;
   } else if (x > 0) {
      z = x;
      u = 0;
      f = 1;
   } else if (x == 0) {
      return 0;
   } else if (x >= -1) {
      l = log(1 - x);
      z = x/(x - 1);
      u = -0.5*l*l;
      f = -1;
   } else {
      z = 1/(1 - x);
      u = 0.5*log(z)*log(x*x*z) - PI6;
      f = 1;
   }

   y = 0.008048124718341*z + 0.008288074835108;
   y = y*z - 0.001481786416153;
   y = y*z - 0.000912777413024;
   y = y*z + 0.005047192127203;
   y = y*z + 0.005300972587634;
   y = y*z + 0.004091615355944;
   y = y*z + 0.004815490327461;
   y = y*z + 0.005966509196748;
   y = y*z + 0.006980881130380;
   y = y*z + 0.008260083434161;
   y = y*z + 0.009997129506220;
   y = y*z + 0.012345919431569;
   y = y*z + 0.015625134938703;
   y = y*z + 0.020408155605916;
   y = y*z + 0.027777774308288;
   y = y*z + 0.040000000124677;
   y = y*z + 0.062500000040762;
   y = y*z + 0.111111111110322;
   y = y*z + 0.249999999999859;
   y = y*z + 1;

   return f*y*z + u;
}
