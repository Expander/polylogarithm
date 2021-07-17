// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the MIT License.
// ====================================================================

#include "Cl2.hpp"
#include "Li2.hpp"

namespace polylogarithm {

/**
 * @brief Clausen function \f$\mathrm{Cl}_2(\theta) = \mathrm{Im}(\mathrm{Li}_2(e^{i\theta}))\f$
 * @param x real angle
 * @return \f$\mathrm{Cl}_2(\theta)\f$
 * @author K.S. Kölbig
 * @note Implementation translated by Alexander Voigt from CERNLIB DCLAUS function C326
 *
 * Journal of Computational and Applied Mathematics 64 (1995) 295-297.
 */
double Cl2(double x) noexcept
{
   const double PI = 3.14159265358979324;
   const double PI2 = 2*PI, PIH = PI/2, RPIH = 2/PI;
   const double A[9] = {
      0.0279528319735756613,
      0.0001763088743898116,
      0.0000012662741461157,
      0.0000000117171818134,
      0.0000000001230064129,
      0.0000000000013952729,
      0.0000000000000166908,
      0.0000000000000002076,
      0.0000000000000000027
   };
   const double B[14] = {
       0.639097088857265341,
      -0.054980569301851716,
      -0.000961261945950606,
      -0.000032054686822550,
      -0.000001329461695426,
      -0.000000062093601824,
      -0.000000003129600656,
      -0.000000000166351954,
      -0.000000000009196527,
      -0.000000000000524004,
      -0.000000000000030580,
      -0.000000000000001820,
      -0.000000000000000110,
      -0.000000000000000007,
   };

   double h = 0;
   double v = std::fmod(std::abs(x), PI2);
   double S = x >= 0 ? 1 : -1;

   if (v > PI) {
      const double p0 = 6.28125;
      const double p1 = 0.0019353071795864769253;
      v = (p0 - v) + p1;
      S = -S;
   }

   if (v == 0 || v == PI) {
      h = 0;
   } else if (v < PIH) {
      const double u = RPIH*v;
      h = 2*u*u - 1;
      const double alfa = h + h;
      double b0 = 0, b1 = 0, b2 = 0;
      for (int i = 8; i >= 0; i--) {
         b0 = A[i] + alfa*b1 - b2;
         b2 = b1;
         b1 = b0;
      }
      h = v*(1 - std::log(v) + 0.5*v*v*(b0 - h*b2));
   } else {
      const double u = RPIH*v - 2;
      h = 2*u*u - 1;
      const double alfa = h + h;
      double b0 = 0, b1 = 0, b2 = 0;
      for (int i = 13; i >= 0; i--) {
         b0 = B[i] + alfa*b1 - b2;
         b2 = b1;
         b1 = b0;
      }
      h = (PI - v)*(b0 - h*b2);
   }

   return S*h;
}

/**
 * @brief Clausen function \f$\mathrm{Cl}_2(\theta) = \mathrm{Im}(\mathrm{Li}_2(e^{i\theta}))\f$ with long double precision
 * @param x real angle
 * @return \f$\mathrm{Cl}_2(\theta)\f$
 */
long double Cl2(long double x) noexcept
{
   const long double PI = 3.14159265358979323846264338327950288L;

   x = std::fmod(std::abs(x), 2*PI);
   const long double s = x >= 0 ? 1 : -1;

   const std::complex<long double> ix(0.0L, x);

   return s*std::imag(Li2(std::exp(ix)));
}

} // namespace polylogarithm
