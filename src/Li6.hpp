// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the GNU Lesser General Public
// License (GNU LGPL) version 3.
// ====================================================================

#pragma once
#include <complex>

namespace polylogarithm {

/// Clausen function with n=6
double Cl6(double);

/// complex polylogarithm with n=5
std::complex<double> Li6(const std::complex<double>&);

} // namespace polylogarithm
