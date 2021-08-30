// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the MIT License.
// ====================================================================

#include "Cl2.hpp"
#include <cmath>

namespace polylogarithm {

/**
 * @brief Clausen function \f$\mathrm{Cl}_2(\theta) = \mathrm{Im}(\mathrm{Li}_2(e^{i\theta}))\f$
 * @param x real angle
 * @return \f$\mathrm{Cl}_2(\theta)\f$
 * @author Alexander Voigt
 * @note Implementation as economized Padé approximation.
 */
double Cl2(double x) noexcept
{
   const double PI = 3.14159265358979324;
   const double PI2 = 2*PI, PIH = PI/2, PI28 = PI*PI/8;
   double sgn = 1;

   if (x < 0) {
      x = -x;
      sgn = -1;
   }

   if (x >= PI2) {
      x = std::fmod(x, PI2);
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
         2.7951565822419270e-02, -8.8865360514541522e-04,
         6.8282348222485902e-06, -7.5276232403566808e-09
      };
      const double Q[] = {
         1.0000000000000000e+00, -3.6904397961160525e-02,
         3.7342870576106476e-04, -8.7460760866531179e-07
      };
      const double y = x*x;
      const double z = y - PI28;
      const double z2 = z*z;
      const double p = P[0] + z * P[1] + z2 * (P[2] + z * P[3]);
      const double q = Q[0] + z * Q[1] + z2 * (Q[2] + z * Q[3]);

      h = x*(1 - std::log(x) + y*p/q/2);
   } else {
      const double P[] = {
         6.4005702446195512e-01, -2.0641655351338783e-01,
         2.4175305223497718e-02, -1.2355955287855728e-03,
         2.5649833551291124e-05, -1.4783829128773320e-07
      };
      const double Q[] = {
         1.0000000000000000e+00, -2.5299102015666356e-01,
         2.2148751048467057e-02, -7.8183920462457496e-04,
         9.5432542196310670e-06, -1.8184302880448247e-08
      };
      const double y = PI - x;
      const double z = y*y - PI28;
      const double z2 = z*z;
      const double z4 = z2*z2;
      const double p = P[0] + z * P[1] + z2 * (P[2] + z * P[3]) +
         z4 * (P[4] + z * P[5]);
      const double q = Q[0] + z * Q[1] + z2 * (Q[2] + z * Q[3]) +
         z4 * (Q[4] + z * Q[5]);

      h = y*p/q;
   }

   return sgn*h;
}

/**
 * @brief Clausen function \f$\mathrm{Cl}_2(\theta) = \mathrm{Im}(\mathrm{Li}_2(e^{i\theta}))\f$ with long double precision
 * @param x real angle
 * @return \f$\mathrm{Cl}_2(\theta)\f$
 * @author K.S. Kölbig
 * @note Implementation translated from CERNLIB DCLAUS function C326
 * and extended to long double precision by Alexander Voigt.
 *
 * Journal of Computational and Applied Mathematics 64 (1995) 295-297.
 */
long double Cl2(long double x) noexcept
{
   const long double PI = 3.14159265358979323846264338327950288L;
   const long double PI2 = 2*PI, PIH = PI/2, RPIH = 2/PI;
   const long double A[19] = {
      0.0279528319735756613494585924765551791L,
      0.0001763088743898115653057636473920103L,
      0.0000012662741461156530021975187159184L,
      0.0000000117171818134392379295428166866L,
      0.0000000001230064128833746922855709386L,
      0.0000000000013952728970012911958374309L,
      0.0000000000000166907761628567345146740L,
      0.0000000000000002076091315145432983502L,
      0.0000000000000000026609198306058056092L,
      0.0000000000000000000349249563561378275L,
      0.0000000000000000000004673313082962865L,
      0.0000000000000000000000063542322337428L,
      0.0000000000000000000000000875698871820L,
      0.0000000000000000000000000012208003299L,
      0.0000000000000000000000000000171890569L,
      0.0000000000000000000000000000002441331L,
      0.0000000000000000000000000000000034940L,
      0.0000000000000000000000000000000000503L,
      0.0000000000000000000000000000000000007L
   };
   const double B[30] = {
      0.6390970888572653413071869135953864197L,
     -0.0549805693018517156397035696498958507L,
     -0.0009612619459506064293859076874070709L,
     -0.0000320546868225504765586825318112711L,
     -0.0000013294616954255450141343828695514L,
     -0.0000000620936018243975194590942773212L,
     -0.0000000031296006563911126723262365339L,
     -0.0000000001663519538192669775933926077L,
     -0.0000000000091965272507194254496027281L,
     -0.0000000000005240037738758450093649037L,
     -0.0000000000000305803841873659454134183L,
     -0.0000000000000018196918249487950988000L,
     -0.0000000000000001100398263196261522324L,
     -0.0000000000000000067451775715424687730L,
     -0.0000000000000000004182784651572477035L,
     -0.0000000000000000000261987180876106127L,
     -0.0000000000000000000016553211620337322L,
     -0.0000000000000000000001053943580037580L,
     -0.0000000000000000000000067562822456482L,
     -0.0000000000000000000000004357492971838L,
     -0.0000000000000000000000000282576222954L,
     -0.0000000000000000000000000018415141361L,
     -0.0000000000000000000000000001205472536L,
     -0.0000000000000000000000000000079233873L,
     -0.0000000000000000000000000000005227403L,
     -0.0000000000000000000000000000000346060L,
     -0.0000000000000000000000000000000022982L,
     -0.0000000000000000000000000000000001531L,
     -0.0000000000000000000000000000000000102L,
     -0.0000000000000000000000000000000000007L
   };

   long double sgn = 1;

   if (x < 0) {
      x = -x;
      sgn = -1;
   }

   if (x >= PI2) {
      x = std::fmod(x, PI2);
   }

   if (x > PI) {
      const long double p0 = 6.28125L;
      const long double p1 = 0.0019353071795864769252867665590057683943L;
      x = (p0 - x) + p1;
      sgn = -sgn;
   }

   if (x == 0 || x == PI) {
      return 0;
   }

   long double h = 0;

   if (x < PIH) {
      const long double u = RPIH*x;
      h = 2*u*u - 1;
      const long double alfa = h + h;
      long double b0 = 0, b1 = 0, b2 = 0;
      for (int i = 18; i >= 0; i--) {
         b0 = A[i] + alfa*b1 - b2;
         b2 = b1;
         b1 = b0;
      }
      h = x*(1 - std::log(x) + x*x*(b0 - h*b2)/2);
   } else {
      const long double u = RPIH*x - 2;
      h = 2*u*u - 1;
      const long double alfa = h + h;
      long double b0 = 0, b1 = 0, b2 = 0;
      for (int i = 29; i >= 0; i--) {
         b0 = B[i] + alfa*b1 - b2;
         b2 = b1;
         b1 = b0;
      }
      h = (PI - x)*(b0 - h*b2);
   }

   return sgn*h;
}

} // namespace polylogarithm
