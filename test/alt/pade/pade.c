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

   double h = 0;

   if (x == 0 || x == PI) {
      h = 0;
   } else if (x < PIH) {
      const double P[] = {
         2.7864056054305254e-02, -1.0323970501300232e-03,
        -9.4562987525464908e-04,  1.6780414902965616e-05,
         8.0197688320475017e-06, -2.9715752436065164e-08,
        -1.0087008387998820e-08
      };
      const double Q[] = {
         1.0000000000000000e+00, -4.4991890469884071e-02,
        -3.8778476053616906e-02,  9.5827859998501856e-04,
         4.2564041369186477e-04, -3.5940355499604070e-06,
        -1.1100637935893686e-06
      };

      const double z = x - PI4;
      const double z2 = z*z;
      const double z4 = z2*z2;
      const double p = P[0] + z * P[1] + z2 * (P[2] + z * P[3]) +
         z4 * (P[4] + z * P[5] + z2 * P[6]);
      const double q = Q[0] + z * Q[1] + z2 * (Q[2] + z * Q[3]) +
         z4 * (Q[4] + z * Q[5] + z2 * Q[6]);

      h = x*(1 - log(x) + x*x*p/q/2);
   } else {
      const double P[] = {
         6.6703664252965825e-1,  6.6614199965545781e-1,
         8.8329092618224283e-2, -9.3825233188409063e-2,
        -2.4538738529998358e-2,  4.2689518564684552e-3,
         1.4780169525251205e-3, -6.8060796972421282e-5,
        -3.1050524815530511e-5,  2.6210699893570077e-7,
         1.7376254458016705e-7
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
      const double p = P[0] + z*P[1] + z2*(P[2] + z*P[3]) +
         z4*(P[4] + z*P[5] + z2*(P[6] + z*P[7])) +
         z8*(P[8] + z*P[9] + z2*P[10]);
      const double q = Q[0] + z*Q[1] + z2*(Q[2] + z*Q[3]) +
         z4*(Q[4] + z*Q[5] + z2*(Q[6] + z*Q[7])) +
         z8*(Q[8] + z*Q[9] + z2*Q[10]);

      h = (PI - x)*p/q;
   }

   return sgn*h;
}
