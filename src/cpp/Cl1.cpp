// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the MIT License.
// ====================================================================

#include "Cl1.hpp"
#include <cmath>
#include <limits>

namespace polylogarithm {

/**
 * @brief Clausen function \f$\operatorname{Cl}_1(\theta) = \operatorname{Re}(\operatorname{Li}_1(e^{i\theta}))\f$
 * @param x real angle
 * @return \f$\operatorname{Cl}_1(\theta)\f$
 * @author Alexander Voigt
 */
double Cl1(double x) noexcept
{
   const double PI = 3.14159265358979324;
   const double PI2 = 2*PI;

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
      return std::numeric_limits<double>::infinity();
   }

   return -std::log(2.0*std::sin(0.5*x));
}

} // namespace polylogarithm
