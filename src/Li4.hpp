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

/// Clausen function with n=4 with long double precision
long double Cl4(long double);

/// complex polylogarithm with n=4
std::complex<double> Li4(const std::complex<double>&);

/// complex polylogarithm with n=4 with long double precision
std::complex<long double> Li4(const std::complex<long double>&);

} // namespace polylogarithm
