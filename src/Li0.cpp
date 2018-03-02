// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the GNU Lesser General Public
// License (GNU LGPL) version 3.
// ====================================================================

#include "Li0.hpp"

namespace polylogarithm {

/**
 * @brief Real polylogarithm \f$\mathrm{Li}_0(x) = x/(1-x)\f$
 * @param x real argument
 * @return \f$\mathrm{Li}_0(x)\f$
 */
double Li0(double x)
{
   return x/(1. - x);
}

/**
 * @brief Complex polylogarithm \f$\mathrm{Li}_0(z) = z/(1-z)\f$
 * @param z complex argument
 * @return \f$\mathrm{Li}_0(z)\f$
 */
std::complex<double> Li0(const std::complex<double>& z)
{
   return z/(1. - z);
}

} // namespace polylogarithm
