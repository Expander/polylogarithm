// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the MIT License.
// ====================================================================

#include "Cl2.hpp"
#include <cmath>

namespace polylogarithm {

/**
 * @brief Clausen function \f$\operatorname{Cl}_2(\theta) = \operatorname{Im}(\operatorname{Li}_2(e^{i\theta}))\f$
 * @param x real angle
 * @return \f$\operatorname{Cl}_2(\theta)\f$
 * @author Alexander Voigt
 * @note Implemented as rational function approximation.
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
         1.3888888888888889e-02, -4.3286930203743071e-04,
         3.2779814789973427e-06, -3.6001540369575084e-09
      };
      const double Q[] = {
         1.0000000000000000e+00, -3.6166589746694121e-02,
         3.6015827281202639e-04, -8.3646182842184428e-07
      };
      const double y = x*x;
      const double y2 = y*y;
      const double p = P[0] + y * P[1] + y2 * (P[2] + y * P[3]);
      const double q = Q[0] + y * Q[1] + y2 * (Q[2] + y * Q[3]);

      h = x*(1 - std::log(x) + y*p/q);
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
 * @brief Clausen function \f$\operatorname{Cl}_2(\theta) = \operatorname{Im}(\operatorname{Li}_2(e^{i\theta}))\f$ with long double precision
 * @param x real angle
 * @return \f$\operatorname{Cl}_2(\theta)\f$
 * @author Alexander Voigt
 *
 * Implemented as an rational function approximation with a maximum
 * error of approximately 3.26e-41 (for long double and quadruple
 * precision).
 */
long double Cl2(long double x) noexcept
{
   const long double PI = 3.14159265358979323846264338327950288L;
   const long double PI2 = 2*PI, PIH = PI/2, PI28 = PI*PI/8;
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
      const long double P[] = {
         1.397578291120963523206040867955323306427e-02L,
        -1.352264019891065415863834380176236559873e-03L,
         5.292882735889014643812912154970234565053e-05L,
        -1.073753987723414895538739538414225390109e-06L,
         1.200707648340635046555652744163248062266e-08L,
        -7.252858952718044684640648391665783926850e-11L,
         2.140267450520462605982610727277758328875e-13L,
        -2.396401118613241903411604342093433593468e-16L,
         4.441828690926415235588391389499684215009e-20L
      };
      const long double Q[] = {
         1.000000000000000000000000000000000000000e+00L,
        -1.018694323414614410071369720193716012304e-01L,
         4.248408782245281612900840467370146443889e-03L,
        -9.337710301347963985908084056584187570954e-05L,
         1.159379163193822053103946363960603543601e-06L,
        -8.083352720393357000801675734774176899515e-09L,
         2.949313240431512997069808854213308209519e-11L,
        -4.742700419624204182400715964695278593777e-14L,
         2.158380636740175386190809152807629331877e-17L
      };
      const long double y = x*x;
      const long double z = y - PI28;
      const long double z2 = z*z;
      const long double z4 = z2*z2;
      const long double z8 = z4*z4;
      const long double p = P[0] + z * P[1] + z2 * (P[2] + z * P[3]) +
         z4 * (P[4] + z * P[5] + z2 * (P[6] + z * P[7])) + z8 * P[8];
      const long double q = Q[0] + z * Q[1] + z2 * (Q[2] + z * Q[3]) +
         z4 * (Q[4] + z * Q[5] + z2 * (Q[6] + z * Q[7])) + z8 * Q[8];

      h = x*(1 - std::log(x) + y*p/q);
   } else {
      const long double P[] = {
         6.400570244619551220929428522356830562481e-01L,
        -4.651631624886004423703445967760673575997e-01L,
         1.487130845262105644024901814213146749895e-01L,
        -2.749665174801454303884783494225610407035e-02L,
         3.251522465413666561950482170352156048203e-03L,
        -2.567438381297475310848635518657180974512e-04L,
         1.372076105130164861564020074178493529151e-05L,
        -4.924179093673498700461153483531075799113e-07L,
         1.153267936031337440182387313169828395860e-08L,
        -1.667578310677508029208023423625588832295e-10L,
         1.348437292247918547169070120217729056878e-12L,
        -5.052245092698477071447850656280954693011e-15L,
         5.600638109466570497480415519182233229048e-18L
      };
      const long double Q[] = {
         1.000000000000000000000000000000000000000e+00L,
        -6.572465772185054284667746526549393897676e-01L,
         1.886234634096976582977630140163583172173e-01L,
        -3.103347567899737687117030083178445406132e-02L,
         3.230860399291224478336071920154030050234e-03L,
        -2.216546569335921728108984951507368387512e-04L,
         1.011949972138985643994631167412420906088e-05L,
        -3.033400935206767852937290458763547850384e-07L,
         5.748454611964843057644023468691231929690e-09L,
        -6.408350048413952604351408631173781906861e-11L,
         3.678584366662951864267349037579031493395e-13L,
        -8.240439699357036167611014086997683837396e-16L,
         3.041049046123062158788159773779755292771e-19L
      };
      const long double y = PI - x;
      const long double z = y*y - PI28;
      const long double z2 = z*z;
      const long double z4 = z2*z2;
      const long double z8 = z4*z4;
      const long double p = P[0] + z * P[1] + z2 * (P[2] + z * P[3]) +
         z4 * (P[4] + z * P[5] + z2 * (P[6] + z * P[7])) +
         z8 * (P[8] + z * P[9] + z2 * (P[10] + z * P[11]) + z4 * P[12]);
      const long double q = Q[0] + z * Q[1] + z2 * (Q[2] + z * Q[3]) +
         z4 * (Q[4] + z * Q[5] + z2 * (Q[6] + z * Q[7])) +
         z8 * (Q[8] + z * Q[9] + z2 * (Q[10] + z * Q[11]) + z4 * Q[12]);

      h = y*p/q;
   }

   return sgn*h;
}

} // namespace polylogarithm
