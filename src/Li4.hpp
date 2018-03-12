// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the GNU Lesser General Public
// License (GNU LGPL) version 3.
// ====================================================================

#pragma once
#include <complex>

namespace polylogarithm {

/// Clausen function with n=4
double Cl4(double);

/// complex polylogarithm with n=4
std::complex<double> Li4(const std::complex<double>&);

} // namespace polylogarithm
