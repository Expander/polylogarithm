// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the MIT License.
// ====================================================================

#include "Cl5.hpp"
#include "Li5.hpp"
#include <complex>

namespace polylogarithm {

/**
 * @brief Clausen function \f$\operatorname{Cl}_5(\theta) = \operatorname{Re}(\operatorname{Li}_5(e^{i\theta}))\f$
 * @param x real angle
 * @return \f$\operatorname{Cl}_5(\theta)\f$
 * @author Alexander Voigt
 * @note Implementation as rational function approximation.
 */
double Cl5(double x) noexcept
{
   const double PI = 3.14159265358979324;
   const double PI2 = 2*PI, PIH = PI/2, PI28 = PI*PI/8;
   const double zeta5 = 1.0369277551433699;

   if (x < 0) {
      x = -x;
   }

   if (x >= PI2) {
      x = std::fmod(x, PI2);
   }

   if (x > PI) {
      const double p0 = 6.28125;
      const double p1 = 0.0019353071795864769253;
      x = (p0 - x) + p1;
   }

   if (x == 0) {
      return zeta5;
   }

   double h = 0;

   if (x < PIH) {
      const double P[] = {
         1.0369277551433699e+00, -6.1354800479984468e-01,
         9.4076401395712763e-02, -9.4056155866704436e-04
      };
      const double Q[] = {
         1.0000000000000000e+00, -1.2073698633244778e-02,
         1.3703409625482991e-05, -1.9701280330628469e-09,
         2.1944550184416500e-11
      };
      const double y = x*x;
      const double y2 = y*y;
      const double p = P[0] + y * P[1] + y2 * (P[2] + y * P[3]);
      const double q = Q[0] + y * Q[1] + y2 * (Q[2] + y * Q[3] + y2 * Q[4]);
      h = p/q - 1./24*y2*std::log(x);
   } else {
      const double P[] = {
         -4.5930112735784898e-01, 4.3720705508867954e-01,
         -7.5895226486465095e-02, 5.2244176912488065e-03,
         -1.5677716622013956e-04, 1.6641624171748576e-06
      };
      const double Q[] = {
          1.0000000000000000e+00, -1.2211486825401188e-01,
          3.8940070749313620e-03, -2.2674805547074318e-05,
         -7.4383354448335299e-08, -3.4131758392216437e-10
      };
      const double y = PI - x;
      const double z = y*y - PI28;
      const double z2 = z*z;
      const double z4 = z2*z2;
      const double p = P[0] + z * P[1] + z2 * (P[2] + z * P[3]) +
         z4 * (P[4] + z * P[5]);
      const double q = Q[0] + z * Q[1] + z2 * (Q[2] + z * Q[3]) +
         z4 * (Q[4] + z * Q[5]);
      h = p/q;
   }

   return h;
}

/**
 * @brief Clausen function \f$\operatorname{Cl}_5(\theta) = \operatorname{Re}(\operatorname{Li}_5(e^{i\theta}))\f$ with long double precision
 * @param x real angle
 * @return \f$\operatorname{Cl}_5(\theta)\f$
 */
long double Cl5(long double x) noexcept
{
   return std::real(Li5(std::polar(1.0L, x)));
}

} // namespace polylogarithm
