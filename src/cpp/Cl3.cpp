// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the MIT License.
// ====================================================================

#include "Cl3.hpp"
#include "Li3.hpp"
#include <complex>

namespace polylogarithm {

/**
 * @brief Clausen function \f$\mathrm{Cl}_3(\theta) = \mathrm{Re}(\mathrm{Li}_3(e^{i\theta}))\f$
 * @param x real angle
 * @return \f$\mathrm{Cl}_3(\theta)\f$
 */
double Cl3(double x) noexcept
{
   return std::real(Li3(std::polar(1.0, x)));
}

/**
 * @brief Clausen function \f$\mathrm{Cl}_3(\theta) = \mathrm{Re}(\mathrm{Li}_3(e^{i\theta}))\f$ with long double precision
 * @param x real angle
 * @return \f$\mathrm{Cl}_3(\theta)\f$
 */
long double Cl3(long double x) noexcept
{
   return std::real(Li3(std::polar(1.0L, x)));
}

} // namespace polylogarithm
