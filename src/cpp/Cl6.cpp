// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the MIT License.
// ====================================================================

#include "Cl6.hpp"
#include "Li6.hpp"
#include <complex>

namespace polylogarithm {

/**
 * @brief Clausen function \f$\operatorname{Cl}_6(\theta) = \operatorname{Im}(\operatorname{Li}_6(e^{i\theta}))\f$
 * @param x real angle
 * @return \f$\operatorname{Cl}_6(\theta)\f$
 * @author Alexander Voigt
 * @note Implemented as a rational function approximation.
 */
double Cl6(double x) noexcept
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
         1.0369277551433699e+00, -2.087195444107175e-01,
         2.0652251045312954e-02, -1.383438138256840e-04
      };
      const double Q[] = {
         1.0000000000000000e+00, -8.0784096827362542e-03,
         5.8074568862993102e-06, -5.1960620033050114e-10
      };
      const double y = x*x;
      const double y2 = y*y;
      const double p = P[0] + y * P[1] + y2 * (P[2] + y * P[3]);
      const double q = Q[0] + y * Q[1] + y2 * (Q[2] + y * Q[3]);
      h = x*(p/q - 1./120*y2*std::log(x));
   } else {
      const double P[] = {
         7.9544504578027050e-01, -1.9255025309738589e-01,
         1.5805208288846591e-02, -5.4175380521534706e-04,
         6.7577493541009068e-06
      };
      const double Q[] = {
         1.0000000000000000e+00, -7.0798422394109274e-02,
         7.1744189715634762e-04,  3.9098747334347093e-06,
         3.5669441618295266e-08,  2.5315391843409925e-10
      };
      const double y = PI - x;
      const double z = y*y - PI28;
      const double z2 = z*z;
      const double z4 = z2*z2;
      const double p = P[0] + z * P[1] + z2 * (P[2] + z * P[3]) +
         z4 * P[4];
      const double q = Q[0] + z * Q[1] + z2 * (Q[2] + z * Q[3]) +
         z4 * (Q[4] + z * Q[5]);
      h = y*p/q;
   }

   return sgn*h;
}

/**
 * @brief Clausen function \f$\operatorname{Cl}_6(\theta) = \operatorname{Im}(\operatorname{Li}_6(e^{i\theta}))\f$ with long double precision
 * @param x real angle
 * @return \f$\operatorname{Cl}_6(\theta)\f$
 */
long double Cl6(long double x) noexcept
{
   return std::imag(Li6(std::polar(1.0L, x)));
}

} // namespace polylogarithm
