// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the MIT License.
// ====================================================================

#include "Cl.hpp"
#include "Li.hpp"
#include <complex>

namespace polylogarithm {

namespace {

bool is_even(int64_t n) { return n % 2 == 0; }

} // anonymous namespace

/**
 * @brief Clausen function \f$\operatorname{Cl}_n(\theta)\f$
 * @param n degree of Clausen function
 * @param x real angle
 * @return \f$\operatorname{Cl}_n(\theta)\f$
 */
double Cl(int64_t n, double x)
{
   const std::complex<double> li = Li(n, std::polar(1.0, x));

   if (is_even(n)) {
      return std::imag(li);
   }

   return std::real(li);
}

} // namespace polylogarithm
