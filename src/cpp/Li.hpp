// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the MIT License.
// ====================================================================

#pragma once
#include <complex>
#include <cstdint>

namespace polylogarithm {

/// complex polylogarithm for arbitrary integer n
std::complex<double> Li(int64_t n, const std::complex<double>&) noexcept;

} // namespace polylogarithm
