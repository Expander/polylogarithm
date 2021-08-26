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
          2.7864056054305254e-02, -1.1518353523811193e-03,
         -8.4418800917714430e-04,  1.8094520214219337e-06,
          6.6281448550088233e-06,  2.3116996576597765e-07,
         -7.9449737962038662e-09, -5.0474502938642823e-10
      };
      const double Q[] = {
         1.0000000000000000e+00, -4.9278355513534521e-02,
        -3.5103839201602514e-02,  4.1409673378170632e-04,
         3.6171200002212000e-04,  8.6992511673171293e-06,
        -9.1614231442660733e-07, -4.7647885659951389e-08
      };

      const double z = x - PI4;
      const double z2 = z*z;
      const double z4 = z2*z2;
      const double p = P[0] + z * P[1] + z2 * (P[2] + z * P[3]) +
         z4 * (P[4] + z * P[5] + z2 * (P[6] + z * P[7]));
      const double q = Q[0] + z * Q[1] + z2 * (Q[2] + z * Q[3]) +
         z4 * (Q[4] + z * Q[5] + z2 * (Q[6] + z * Q[7]));

      h = x*(1 - log(x) + x*x*p/q/2);
   } else {
      const double P[] = {
         6.6703664252965825e-01,  1.2253096312594806e+00,
         7.2049312875025059e-01,  3.1068589986833602e-02,
        -1.0924727599310961e-01, -2.5092151007588517e-02,
         4.9306731929865383e-03,  1.7068884136587118e-03,
        -7.2196440547828341e-05, -3.8570269835412217e-05,
         2.0126671630631965e-07,  2.3086886135348935e-07
      };
      const double Q[] = {
         1.0000000000000000e+00,  1.7356582515535528e+00,
         9.7304146556104840e-01,  6.1554229775992906e-02,
        -1.1068913603979054e-01, -2.4089856010563872e-02,
         3.5143806119270642e-03,  1.1044099641843965e-03,
        -3.1518129467240003e-05, -1.4600588443999088e-05,
         3.3309035195403625e-08,  2.8909514867371887e-08
      };
      const double z = x - PI34;
      const double z2 = z*z;
      const double z4 = z2*z2;
      const double z8 = z4*z4;
      const double p = P[0] + z*P[1] + z2*(P[2] + z*P[3]) +
         z4*(P[4] + z*P[5] + z2*(P[6] + z*P[7])) +
         z8*(P[8] + z*P[9] + z2*(P[10] + z*P[11]));
      const double q = Q[0] + z*Q[1] + z2*(Q[2] + z*Q[3]) +
         z4*(Q[4] + z*Q[5] + z2*(Q[6] + z*Q[7])) +
         z8*(Q[8] + z*Q[9] + z2*(Q[10] + z*Q[11]));

      h = (PI - x)*p/q;
   }

   return sgn*h;
}
