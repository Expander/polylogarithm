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

/// real dilogarithm
double Li2(double);

/// complex dilogarithm
std::complex<double> Li2(const std::complex<double>&);

/// real trilogarithm
double Li3(double);

// Clausen function
double Cl2(double);

} // namespace dilogarithm

#endif
