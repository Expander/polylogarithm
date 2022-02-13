// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the MIT License.
// ====================================================================

#include "harmonic.hpp"

namespace polylogarithm {

namespace {

} // anonymous namespace

/// harmonic number n
double harmonic(int64_t n) noexcept
{
   double sum = 0;

   for (int64_t h = 1; h <= n; h++) {
      sum += 1./h;
   }

   return sum;
}

} // namespace polylogarithm
