// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the MIT License.
// ====================================================================

#include "Cl6.hpp"
#include "Li6.hpp"
#include <cmath>

namespace polylogarithm {

/**
 * @brief Clausen function \f$\mathrm{Cl}_6(\theta) = \mathrm{Im}(\mathrm{Li}_6(e^{i\theta}))\f$
 * @param x real angle
 * @return \f$\mathrm{Cl}_6(\theta)\f$
 */
double Cl6(double x) noexcept
{
   const std::complex<double> ix(0.0, x);
   return std::imag(Li6(std::exp(ix)));
}

/**
 * @brief Clausen function \f$\mathrm{Cl}_6(\theta) = \mathrm{Im}(\mathrm{Li}_6(e^{i\theta}))\f$ with long double precision
 * @param x real angle
 * @return \f$\mathrm{Cl}_6(\theta)\f$
 */
long double Cl6(long double x) noexcept
{
   const std::complex<long double> ix(0.0L, x);
   return std::imag(Li6(std::exp(ix)));
}

} // namespace polylogarithm
