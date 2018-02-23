// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the GNU Lesser General Public
// License (GNU LGPL) version 3.
// ====================================================================

#ifndef POLYLOGARITHM_H
#define POLYLOGARITHM_H

#include <complex>

namespace polylogarithm {

namespace {
   const int version_major = 2;
   const int version_minor = 0;
   const int version_patch = 0;
}

/// real polylogarithm with n=0
double Li0(double);

/// complex polylogarithm with n=0
std::complex<double> Li0(const std::complex<double>&);

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

/// complex polylogarithm with n=4
std::complex<double> Li4(const std::complex<double>&);

/// complex polylogarithm with n=5
std::complex<double> Li5(const std::complex<double>&);

/// complex polylogarithm with n=5
std::complex<double> Li6(const std::complex<double>&);

// Clausen function
double Cl2(double);

} // namespace polylogarithm

#endif
