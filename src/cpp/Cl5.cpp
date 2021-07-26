// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the MIT License.
// ====================================================================

#include "Cl5.hpp"
#include "Li5.hpp"
#include <complex>

namespace polylogarithm {

/**
 * @brief Clausen function \f$\mathrm{Cl}_5(\theta) = \mathrm{Re}(\mathrm{Li}_5(e^{i\theta}))\f$
 * @param x real angle
 * @return \f$\mathrm{Cl}_5(\theta)\f$
 */
double Cl5(double x) noexcept
{
   return std::real(Li5(std::polar(1.0, x)));
}

/**
 * @brief Clausen function \f$\mathrm{Cl}_5(\theta) = \mathrm{Re}(\mathrm{Li}_5(e^{i\theta}))\f$ with long double precision
 * @param x real angle
 * @return \f$\mathrm{Cl}_5(\theta)\f$
 */
long double Cl5(long double x) noexcept
{
   return std::real(Li5(std::polar(1.0L, x)));
}

} // namespace polylogarithm
