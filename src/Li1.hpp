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

/// Clausen function with n=1 with long double precision
long double Cl1(long double);

/// real polylogarithm with n=1
double Li1(double);

/// real polylogarithm with n=1 with long double precision
long double Li1(long double);

/// complex polylogarithm with n=1
std::complex<double> Li1(const std::complex<double>&);

/// complex polylogarithm with n=1 with long double precision
std::complex<long double> Li1(const std::complex<long double>&);

} // namespace polylogarithm
