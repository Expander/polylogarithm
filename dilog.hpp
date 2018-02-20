// ====================================================================
// This file is part of Dilogarithm.
//
// Dilogarithm is licenced under the GNU Lesser General Public License
// (GNU LGPL) version 3.
// ====================================================================

#ifndef DILOG_H
#define DILOG_H

#include <complex>

namespace dilogarithm {

/// real polylogarithm with n=1
double Li1(double);

/// complex polylogarithm with n=1
std::complex<double> Li1(const std::complex<double>&);

/// real dilogarithm
double Li2(double);

/// complex dilogarithm
std::complex<double> Li2(const std::complex<double>&);

/// complex trilogarithm
std::complex<double> Li3(const std::complex<double>&);

// Clausen function
double Cl2(double);

} // namespace dilogarithm

#endif
