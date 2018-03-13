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
   std::complex<double> clog(std::complex<double> z) noexcept {
      if (std::real(z) == 0.0) z.real(0.0);
      if (std::imag(z) == 0.0) z.imag(0.0);
      return std::log(z);
   }
} // anonymous namespace

/**
 * @brief Clausen function \f$\mathrm{Cl}_1(\theta) = \mathrm{Re}(\mathrm{Li}_1(e^{i\theta}))\f$
 * @param x real angle
 * @return \f$\mathrm{Cl}_1(\theta)\f$
 */
double Cl1(double x)
{
   using std::exp;
   const std::complex<double> i(0.,1.);

   return std::real(Li1(exp(i*x)));
}

/**
 * @brief Real polylogarithm \f$\mathrm{Li}_1(x) = -\log(1-x)\f$
 * @param x real argument
 * @return \f$\mathrm{Li}_1(x)\f$
 */
double Li1(double x)
{
   return -std::log(1. - x);
}

/**
 * @brief Complex polylogarithm \f$\mathrm{Li}_1(z) = -\log(1-z)\f$
 * @param z complex argument
 * @return \f$\mathrm{Li}_1(z)\f$
 */
std::complex<double> Li1(const std::complex<double>& z)
{
   return -clog(1. - z);
}

} // namespace polylogarithm
