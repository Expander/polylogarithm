#include <math.h>

/**
 * @brief Clausen function \f$\mathrm{Cl}_2(\theta) = \mathrm{Im}(\mathrm{Li}_2(e^{i\theta}))\f$
 * @param x real angle
 * @return \f$\mathrm{Cl}_2(\theta)\f$
 * @author Alexander Voigt
 *
 * Implementation as series expansion from Jiming Wu, Xiaoping Zhang,
 * Dongjie Liu: "An efficient calculation of the Clausen functions
 * Cl_n(0)(n >= 2)"
 */
double clausen_2_wu(double x)
{
   const double PI = 3.14159265358979324;
   const double PI2 = 2*PI;
   const double C[12] = {
      1.0,
      -0.027777777777777778,
      -0.00027777777777777778,
      -4.7241118669690098e-06,
      -9.1857730746619636e-08,
      -1.8978869988970999e-09,
      -4.0647616451442255e-11,
      -8.9216910204564526e-13,
      -1.9939295860721076e-14,
      -4.5189800296199182e-16,
      -1.0356517612181247e-17,
      -2.3952186210261867e-19
   };

   double sgn = 1;

   if (x < 0) {
      x = -x;
      sgn = -1;
   }

   if (x >= PI2) {
      x = fmod(x, PI2);
   }

   if (x > PI) {
      const double p0 = 6.28125;
      const double p1 = 0.0019353071795864769253;
      x = (p0 - x) + p1;
      sgn = -sgn;
   }

   if (x == 0 || x == PI) {
      return 0;
   }

   const double x2 = x*x;
   const double sum =
      - log(2*sin(x/2))
      + C[0] +
      + x2*(C[1]
      + x2*(C[2]
      + x2*(C[3]
      + x2*(C[4]
      + x2*(C[5]
      + x2*(C[6]
      + x2*(C[7]
      + x2*(C[8]
      + x2*(C[9]
      + x2*(C[10]
      + x2*(C[11])))))))))));

   return sgn*x*sum;
}
