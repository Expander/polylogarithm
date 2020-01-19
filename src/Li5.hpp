// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the GNU Lesser General Public
// License (GNU LGPL) version 3.
// ====================================================================

#pragma once
#include <complex>

namespace polylogarithm {

/// Clausen function with n=5
double Cl5(double);

/// Clausen function with n=5 with long double precision
long double Cl5(long double);

/// complex polylogarithm with n=5
std::complex<double> Li5(const std::complex<double>&);

/// complex polylogarithm with n=5 with long double precision
std::complex<long double> Li5(const std::complex<long double>&);

} // namespace polylogarithm
