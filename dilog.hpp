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
double dilog(double);

/// complex dilogarithm
std::complex<double> dilog(const std::complex<double>&);

// Clausen function
double clausen_2(double);

} // namespace dilogarithm

#endif
