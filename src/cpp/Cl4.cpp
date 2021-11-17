// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the MIT License.
// ====================================================================

#include "Cl4.hpp"
#include "Li4.hpp"
#include <cmath>
#include <complex>

namespace polylogarithm {

/**
 * @brief Clausen function \f$\mathrm{Cl}_4(\theta) = \mathrm{Im}(\mathrm{Li}_4(e^{i\theta}))\f$
 * @param x real angle
 * @return \f$\mathrm{Cl}_4(\theta)\f$
 * @author Alexander Voigt
 * @note Implemented as economized Padé approximation.
 */
double Cl4(double x) noexcept
{
   const double PI = 3.14159265358979324;
   const double PI2 = 2*PI, PIH = PI/2, PI28 = PI*PI/8;
   const double zeta3 = 1.2020569031595943;
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
         -3.0641482939025622e-01,  6.1700119337868541e-03,
         -2.0244370294666391e-05, -3.1997168417491939e-08
      };
      const double Q[] = {
         1.0000000000000000e+00, -2.2415973860234228e-02,
         1.1164184654978598e-04, -6.3541742491831717e-08
      };
      const double y = x*x;
      const double z = y - PI28;
      const double z2 = z*z;
      const double p = P[0] + z * P[1] + z2 * (P[2] + z * P[3]);
      const double q = Q[0] + z * Q[1] + z2 * (Q[2] + z * Q[3]);

      h = x*(zeta3 + y*(p/q + std::log(x)/6));
   } else {
      const double P[] = {
         7.6223911686491336e-01, -2.4339587368267260e-01,
         2.8715364937979943e-02, -1.5368612510964667e-03,
         3.6261044225761673e-05, -2.8557977333851308e-07
      };
      const double Q[] = {
         1.0000000000000000e+00, -1.7465715261403233e-01,
         9.5439417991615653e-03, -1.7325070821666274e-04,
         5.9283675098376635e-07,  9.4127575773361230e-10
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
 * @brief Clausen function \f$\mathrm{Cl}_4(\theta) = \mathrm{Im}(\mathrm{Li}_4(e^{i\theta}))\f$ with long double precision
 * @param x real angle
 * @return \f$\mathrm{Cl}_4(\theta)\f$
 */
long double Cl4(long double x) noexcept
{
   return std::imag(Li4(std::polar(1.0L, x)));
}

} // namespace polylogarithm
