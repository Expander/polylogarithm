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
   const double eps_d = 10.0*std::numeric_limits<double>::epsilon();
   const long double eps_ld = 10.0L*std::numeric_limits<long double>::epsilon();
   const double inf_d = std::numeric_limits<double>::infinity();
   const long double inf_ld = std::numeric_limits<long double>::infinity();

   template <typename T>
   bool is_close(T a, T b, T eps)
   {
      return std::abs(a - b) < eps;
   }

   template <typename T>
   bool is_close(const std::complex<T>& a, T b, T eps)
   {
      return std::abs(std::real(a) - b) < eps && std::abs(std::imag(a)) < eps;
   }
} // anonymous namespace

/**
 * @brief Real polylogarithm \f$\mathrm{Li}_0(x) = x/(1-x)\f$
 * @param x real argument
 * @return \f$\mathrm{Li}_0(x)\f$
 */
double Li0(double x)
{
   if (is_close(x, 1.0, eps_d)) {
      return inf_d;
   }

   return x/(1.0 - x);
}

/**
 * @brief Real polylogarithm \f$\mathrm{Li}_0(x) = x/(1-x)\f$ with long double precision
 * @param x real argument
 * @return \f$\mathrm{Li}_0(x)\f$
 */
long double Li0(long double x)
{
   if (is_close(x, 1.0L, eps_ld)) {
      return inf_ld;
   }

   return x/(1.0L - x);
}

/**
 * @brief Complex polylogarithm \f$\mathrm{Li}_0(z) = z/(1-z)\f$
 * @param z complex argument
 * @return \f$\mathrm{Li}_0(z)\f$
 */
std::complex<double> Li0(const std::complex<double>& z)
{
   if (is_close(z, 1.0, eps_d)) {
      return {inf_d, inf_d};
   }

   return z/(1.0 - z);
}

/**
 * @brief Complex polylogarithm \f$\mathrm{Li}_0(z) = z/(1-z)\f$ with long double precision
 * @param z complex argument
 * @return \f$\mathrm{Li}_0(z)\f$
 */
std::complex<long double> Li0(const std::complex<long double>& z)
{
   if (is_close(z, 1.0L, eps_ld)) {
      return {inf_ld, inf_ld};
   }

   return z/(1.0L - z);
}

} // namespace polylogarithm
