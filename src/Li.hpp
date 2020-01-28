// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the GNU Lesser General Public
// License (GNU LGPL) version 3.
// ====================================================================

#pragma once
#include <complex>

namespace polylogarithm {

/// Clausen function for arbitrary integer n
double Cl(int64_t, double);

/// complex polylogarithm for arbitrary integer n
std::complex<double> Li(int64_t n, const std::complex<double>&);

} // namespace polylogarithm
