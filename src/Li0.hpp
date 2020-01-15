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

/// real polylogarithm with n=0 with long double precision
long double Li0(long double);

/// complex polylogarithm with n=0
std::complex<double> Li0(const std::complex<double>&);

/// complex polylogarithm with n=0 with long double precision
std::complex<long double> Li0(const std::complex<long double>&);

} // namespace polylogarithm
