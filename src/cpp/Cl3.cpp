// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the MIT License.
// ====================================================================

#include "Cl3.hpp"
#include "Li3.hpp"
#include <cmath>
#include <complex>

namespace polylogarithm {

/**
 * @brief Clausen function \f$\mathrm{Cl}_3(\theta) = \mathrm{Re}(\mathrm{Li}_3(e^{i\theta}))\f$
 * @param x real angle
 * @return \f$\mathrm{Cl}_3(\theta)\f$
 * @author Alexander Voigt
 * @note Implementation as economized Padé approximation.
 */
double Cl3(double x) noexcept
{
   const double PI = 3.14159265358979324;
   const double PI2 = 2*PI, PIH = PI/2, PI28 = PI*PI/8;
   const double zeta3 = 1.2020569031595943;

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
      return zeta3;
   }

   double h = 0;

   if (x < PIH) {
      const double P[] = {
         -7.5430148591242361e-01,  1.6121940167854339e-02,
         -3.7484056212140535e-05, -2.5191292110169198e-07
      };
      const double Q[] = {
         1.0000000000000000e+00, -2.6015033560727570e-02,
         1.5460630299236049e-04, -1.0987530650923219e-07
      };
      const double y = x*x;
      const double z = y - PI28;
      const double z2 = z*z;
      const double p = P[0] + z * P[1] + z2 * (P[2] + z * P[3]);
      const double q = Q[0] + z * Q[1] + z2 * (Q[2] + z * Q[3]);
      h = zeta3 + y*(p/q + std::log(x)/2);
   } else {
      const double P[] = {
         -4.9017024647634973e-01, 4.1559155224660940e-01,
         -7.9425531417806701e-02, 5.9420152260602943e-03,
         -1.8302227163540190e-04, 1.8027408929418533e-06
      };
      const double Q[] = {
         1.0000000000000000e+00, -1.9495887541644712e-01,
         1.2059410236484074e-02, -2.5235889467301620e-04,
         1.0199322763377861e-06,  1.9612106499469264e-09
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
 * @brief Clausen function \f$\mathrm{Cl}_3(\theta) = \mathrm{Re}(\mathrm{Li}_3(e^{i\theta}))\f$ with long double precision
 * @param x real angle
 * @return \f$\mathrm{Cl}_3(\theta)\f$
 */
long double Cl3(long double x) noexcept
{
   return std::real(Li3(std::polar(1.0L, x)));
}

} // namespace polylogarithm
