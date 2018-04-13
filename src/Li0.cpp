// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the GNU Lesser General Public
// License (GNU LGPL) version 3.
// ====================================================================

#include "Li0.hpp"
#include <cmath>
#include <limits>

namespace polylogarithm {

namespace {
   const double epsilon = std::pow(10., -std::floor(std::numeric_limits<double>::digits10));
   const double inf = std::numeric_limits<double>::infinity();

   bool is_close(double a, double b, double eps = epsilon)
   {
      return std::abs(a - b) < eps;
   }

   bool is_close(const std::complex<double>& a, const std::complex<double>& b,
                 double eps = epsilon)
   {
      return std::abs(std::real(a) - std::real(b)) < eps &&
             std::abs(std::imag(a) - std::imag(b)) < eps;
   }
} // anonymous namespace

/**
 * @brief Real polylogarithm \f$\mathrm{Li}_0(x) = x/(1-x)\f$
 * @param x real argument
 * @return \f$\mathrm{Li}_0(x)\f$
 */
double Li0(double x)
{
   if (is_close(x, 1.))
      return inf;

   return x/(1. - x);
}

/**
 * @brief Complex polylogarithm \f$\mathrm{Li}_0(z) = z/(1-z)\f$
 * @param z complex argument
 * @return \f$\mathrm{Li}_0(z)\f$
 */
std::complex<double> Li0(const std::complex<double>& z)
{
   if (is_close(z, {1.,0.}))
      return {inf, inf};

   return z/(1. - z);
}

} // namespace polylogarithm
