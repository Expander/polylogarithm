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
   const double PI2 = 2*PI, PIH = PI/2, RPIH = 2/PI;
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
          0.027864056054305254, -0.0011518353523811193,
         -0.0008441880091771443, 1.8094520214219337e-6,
          6.6281448550088233e-6, 2.3116996576597765e-7,
         -7.9449737962038662e-9, -5.0474502938642823e-10
      };
      const double Q[] = {
         1., -0.049278355513534521,
         -0.035103839201602514, 0.00041409673378170632,
         0.00036171200002212, 8.6992511673171293e-6,
         -9.1614231442660733e-7, -4.7647885659951389e-8};

      const double z = x - PIH/2;
      const double z2 = z*z;
      const double z4 = z2*z2;
      const double p = P[0] + z * P[1] + z2 * (P[2] + z * P[3]) +
         z4 * (P[4] + z * P[5] + z2 * (P[6] + z * P[7]));
      const double q = Q[0] + z * Q[1] + z2 * (Q[2] + z * Q[3]) +
         z4 * (Q[4] + z * Q[5] + z2 * (Q[6] + z * Q[7]));

      h = x*(1 - log(x) + x*x*p/q/2);
   } else {
      /* @todo: improve precision */
      const double P[] = {
         0.66703664252965849, 0.85935948413199535,
         0.52264499054693664, 0.16483396470454956,
         0.0011202318460937339, -0.021627581091508326,
         -0.0076882879112172301, -0.00011906461695951715,
         0.00051607758453434685, 0.000068216104614467083,
         -0.000016281772164610036, -1.6848812010005932e-6,
         2.6890142240624932e-7
      };
      const double Q[] = {
         1., 1.187037495518838,
         0.73200188213651604, 0.24881464958122803,
         0.02237477346709609, -0.019440158049282673,
         -0.0081044741661252977, -0.00042872698131159565,
         0.00038998159876750839, 0.000044928990350554166,
         -9.7986918698359905e-6, -3.7158127315397317e-7,
         5.8923843438724757e-8
      };
      const double z = x - 3*PIH/2;
      const double z2 = z*z;
      const double z4 = z2*z2;
      const double z8 = z4*z4;
      const double p = P[0] + z*P[1] + z2*(P[2] + z*P[3]) +
         z4*(P[4] + z*P[5] + z2*(P[6] + z*P[7])) +
         z8*(P[8] + z*P[9] + z2*(P[10] + z*P[11]) + z4*P[12]);
      const double q = Q[0] + z*Q[1] + z2*(Q[2] + z*Q[3]) +
         z4*(Q[4] + z*Q[5] + z2*(Q[6] + z*Q[7])) +
         z8*(Q[8] + z*Q[9] + z2*(Q[10] + z*Q[11]) + z4*Q[12]);

      h = (PI - x)*p/q;
   }

   return sgn*h;
}
