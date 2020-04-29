// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the GNU Lesser General Public
// License (GNU LGPL) version 3.
// ====================================================================

#include "Cl1.hpp"
#include "Li1.hpp"
#include <cmath>

namespace polylogarithm {

/**
 * @brief Clausen function \f$\mathrm{Cl}_1(\theta) = \mathrm{Re}(\mathrm{Li}_1(e^{i\theta}))\f$
 * @param x real angle
 * @return \f$\mathrm{Cl}_1(\theta)\f$
 */
double Cl1(double x)
{
   const std::complex<double> i(0.0, 1.0);

   return std::real(Li1(std::exp(i*x)));
}

/**
 * @brief Clausen function \f$\mathrm{Cl}_1(\theta) = \mathrm{Re}(\mathrm{Li}_1(e^{i\theta}))\f$ with long double precision
 * @param x real angle
 * @return \f$\mathrm{Cl}_1(\theta)\f$
 */
long double Cl1(long double x)
{
   const std::complex<long double> i(0.0L, 1.0L);

   return std::real(Li1(std::exp(i*x)));
}

} // namespace polylogarithm
