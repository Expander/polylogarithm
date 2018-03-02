// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the GNU Lesser General Public
// License (GNU LGPL) version 3.
// ====================================================================

#pragma once
#include <complex>

namespace polylogarithm {

/// real polylogarithm with n=0
double Li0(double);

/// complex polylogarithm with n=0
std::complex<double> Li0(const std::complex<double>&);

} // namespace polylogarithm
