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

      const double z = x*x - PI4;
      const double z2 = z*z;
      const double p = P[0] + z * P[1] + z2 * (P[2] + z * P[3]);
      const double q = Q[0] + z * Q[1] + z2 * (Q[2] + z * Q[3]);

      h = x*(1 - log(x) + x*x*p/q/2);
   } else {
      const double P[] = {
         6.6703664252965825e-01,  6.6614199965545781e-01,
         8.8329092618224283e-02, -9.3825233188409063e-02,
        -2.4538738529998358e-02,  4.2689518564684552e-03,
         1.4780169525251205e-03, -6.8060796972421282e-05,
        -3.1050524815530511e-05,  2.6210699893570077e-07,
         1.7376254458016705e-07
      };
      const double Q[] = {
         1.0000000000000000e+00,  8.9737202265397515e-01,
         1.1022865557546485e-01, -9.5881689003093915e-02,
        -2.2242106921505683e-02,  3.0797298590089902e-03,
         9.2028848593374794e-04, -2.9945345568046459e-05,
        -1.1356385746761059e-05,  4.1088867204827090e-08,
         2.1073459001575500e-08
      };
      const double z = x - PI34;
      const double z2 = z*z;
      const double z4 = z2*z2;
      const double z8 = z4*z4;
      const double p = P[0] + z * P[1] + z2 * (P[2] + z * P[3]) +
         z4 * (P[4] + z * P[5] + z2 * (P[6] + z * P[7])) +
         z8 * (P[8] + z * P[9] + z2 * P[10]);
      const double q = Q[0] + z * Q[1] + z2 * (Q[2] + z * Q[3]) +
         z4 * (Q[4] + z * Q[5] + z2 * (Q[6] + z * Q[7])) +
         z8 * (Q[8] + z * Q[9] + z2 * Q[10]);

      h = (PI - x)*p/q;
   }

   return sgn*h;
}
