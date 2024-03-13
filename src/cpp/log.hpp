// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the MIT License.
// ====================================================================

#pragma once

#include <cmath>
#include <complex>

namespace polylogarithm {


template <typename T>
std::complex<T> log1p(const std::complex<T>& z) noexcept
{
   const std::complex<T> u = T(1) + z;

   if (std::real(u) == T(1) && std::imag(u) == T(0)) {
      return z;
   } else if (std::real(u) <= T(0)) {
      return log(u);
   }

   return std::log(u)*(z/(u - T(1)));
}


} // namespace polylogarithm
