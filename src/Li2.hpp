// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the GNU Lesser General Public
// License (GNU LGPL) version 3.
// ====================================================================

#pragma once
#include <complex>

namespace polylogarithm {

/// Clausen function with n=2
double Cl2(double);

/// Clausen function with n=2 with long double precision
long double Cl2(long double);

/// real polylogarithm with n=2 (dilogarithm)
double Li2(double);

/// real polylogarithm with n=2 (dilogarithm) with long double precision
long double Li2(long double);

/// complex polylogarithm with n=2 (dilogarithm)
std::complex<double> Li2(const std::complex<double>&);

/// complex polylogarithm with n=2 (dilogarithm) with long double precision
std::complex<long double> Li2(const std::complex<long double>&);

} // namespace polylogarithm
