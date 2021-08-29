#include <math.h>

/**
 * @brief Clausen function \f$\mathrm{Cl}_2(\theta) = \mathrm{Im}(\mathrm{Li}_2(e^{i\theta}))\f$
 * @param x real angle
 * @return \f$\mathrm{Cl}_2(\theta)\f$
 * @author Alexander Voigt
 * @note Implementation as economized Padé approximation.
 */
double clausen_2_pade(double x)
{
   const double PI = 3.14159265358979324;
   const double PI2 = 2*PI, PIH = PI/2, PI4 = PI/4, PI34 = 3*PI/4;
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

   double h = 0;

   if (x < PIH) {
      const double P[] = {
         2.7887843484730669e-02, -8.7658005565574189e-04,
         6.6596845091611126e-06, -7.2581927220626547e-09
      };
      const double Q[] = {
         1.0000000000000000e+00, -3.6502848628035564e-02,
         3.6543318130278056e-04, -8.4709594747969947e-07
      };
      const double y = x*x;
      const double z = y - PI4;
      const double z2 = z*z;
      const double p = P[0] + z * P[1] + z2 * (P[2] + z * P[3]);
      const double q = Q[0] + z * Q[1] + z2 * (Q[2] + z * Q[3]);

      h = x*(1 - log(x) + y*p/q/2);
   } else {

      const double P[] = {
         5.8843335892981230e-01, -2.5713541030583132e-01,
         4.3178185208549071e-02, -3.4747514300508753e-03,
         1.3555812078275216e-04, -2.2560217492329007e-06,
         1.0674776164945825e-08
      };
      const double Q[] = {
         1.0000000000000000e+00, -3.5610339469990070e-01,
         4.7166170516951617e-02, -2.8540514457813485e-03,
         7.7361534856262963e-05, -7.6064119916616840e-07,
         1.1922724268534526e-09
      };
      const double y = PI - x;
      const double z = y*y - PI34;
      const double z2 = z*z;
      const double z4 = z2*z2;
      const double p = P[0] + z * P[1] + z2 * (P[2] + z * P[3]) +
         z4 * (P[4] + z * P[5] + z2 * (P[6]));
      const double q = Q[0] + z * Q[1] + z2 * (Q[2] + z * Q[3]) +
         z4 * (Q[4] + z * Q[5] + z2 * (Q[6]));

      h = y*p/q;
   }

   return sgn*h;
}
