// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the MIT License.
// ====================================================================

#pragma once
#include <complex>

namespace polylogarithm {

/// complex polylogarithm with n=4
std::complex<double> Li4(const std::complex<double>&) noexcept;

/// complex polylogarithm with n=4 with long double precision
std::complex<long double> Li4(const std::complex<long double>&) noexcept;

} // namespace polylogarithm
