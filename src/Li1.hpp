// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the GNU Lesser General Public
// License (GNU LGPL) version 3.
// ====================================================================

#pragma once
#include <complex>

namespace polylogarithm {

/// Clausen function with n=1
double Cl1(double);

/// real polylogarithm with n=1
double Li1(double);

/// complex polylogarithm with n=1
std::complex<double> Li1(const std::complex<double>&);

} // namespace polylogarithm
