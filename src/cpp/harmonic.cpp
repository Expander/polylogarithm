// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the MIT License.
// ====================================================================

#include "harmonic.hpp"
#include <cmath>
#include <limits>

namespace polylogarithm {

namespace {

/// digamma for integer n > 0, following
/// [K.S. Kölbig: Programs for computing the logarithm of the gamma
/// function, and the digamma function, for complex argument, Computer
/// Physics Communications, Volume 4, Issue 2, 1972, Pages 221-226, ISSN
/// 0010-4655, https://doi.org/10.1016/0010-4655(72)90012-4]
double digamma(int64_t n) noexcept
{
   // Table[BernoulliB[2n]/(2 n), {n,1,8}]
   const double c[] = {
      0.083333333333333333, -0.0083333333333333333,  0.0039682539682539683,
     -0.0041666666666666667, 0.0075757575757575758, -0.021092796092796093,
      0.083333333333333333, -0.44325980392156863
   };

   if (n <= 0) {
      return std::numeric_limits<double>::quiet_NaN();
   }

   double res = 0;

   // map potentially small n to n >= 7
   if (n < 7) { // recurrence formula
      for (int64_t nu = 1; nu < 7 - n; ++nu) {
         res -= 1.0/(n + nu);
      }
      res -= 1.0/n;
      n = 7.0;
   }

   const double t = 1.0/n;
   const double t2 = t*t;

   return res + std::log(n) - 0.5*t
      - t2*(c[0] + t2*(c[1] + t2*(c[2] + t2*(c[3] + t2*(c[4] + t2*(c[5] + t2*(c[6] + t2*c[7])))))));
}

} // anonymous namespace

/// harmonic number n
double harmonic(int64_t n) noexcept
{
   if (n <= 0) {
      return std::numeric_limits<double>::quiet_NaN();
   } else if (n < 20) {
      double sum = 1;
      for (int64_t k = 2; k <= n; ++k) {
         sum += 1.0/k;
      }
      return sum;
   } else {
      const double eulergamma = 0.57721566490153286;
      return eulergamma + digamma(n + 1);
   }
}

} // namespace polylogarithm
