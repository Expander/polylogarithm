// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the MIT License.
// ====================================================================

#include "Cl4.hpp"
#include "Li4.hpp"
#include <complex>

namespace polylogarithm {

/**
 * @brief Clausen function \f$\mathrm{Cl}_4(\theta) = \mathrm{Im}(\mathrm{Li}_4(e^{i\theta}))\f$
 * @param x real angle
 * @return \f$\mathrm{Cl}_4(\theta)\f$
 */
double Cl4(double x) noexcept
{
   return std::imag(Li4(std::polar(1.0, x)));
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
