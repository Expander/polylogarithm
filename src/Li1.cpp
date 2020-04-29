// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the GNU Lesser General Public
// License (GNU LGPL) version 3.
// ====================================================================

#include "Li1.hpp"
#include <cmath>

namespace polylogarithm {

namespace {
   // converts -0.0 to 0.0
   template <typename T>
   std::complex<T> clog(std::complex<T> z) noexcept {
      if (std::real(z) == T(0)) { z.real(T(0)); }
      if (std::imag(z) == T(0)) { z.imag(T(0)); }
      return std::log(z);
   }
} // anonymous namespace

/**
 * @brief Real polylogarithm \f$\mathrm{Li}_1(x) = -\log(1-x)\f$
 * @param x real argument
 * @return \f$\mathrm{Li}_1(x)\f$
 */
double Li1(double x)
{
   return -std::log(1.0 - x);
}

/**
 * @brief Real polylogarithm \f$\mathrm{Li}_1(x) = -\log(1-x)\f$ with long double precision
 * @param x real argument
 * @return \f$\mathrm{Li}_1(x)\f$
 */
long double Li1(long double x)
{
   return -std::log(1.0L - x);
}

/**
 * @brief Complex polylogarithm \f$\mathrm{Li}_1(z) = -\log(1-z)\f$
 * @param z complex argument
 * @return \f$\mathrm{Li}_1(z)\f$
 */
std::complex<double> Li1(const std::complex<double>& z)
{
   return -clog(1.0 - z);
}

/**
 * @brief Complex polylogarithm \f$\mathrm{Li}_1(z) = -\log(1-z)\f$ with long double precision
 * @param z complex argument
 * @return \f$\mathrm{Li}_1(z)\f$
 */
std::complex<long double> Li1(const std::complex<long double>& z)
{
   return -clog(1.0L - z);
}

} // namespace polylogarithm
