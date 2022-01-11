// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the MIT License.
// ====================================================================

#include "Cl.hpp"
#include "Cl1.hpp"
#include "Li.hpp"
#include <complex>
#include <limits>

namespace polylogarithm {

namespace {

constexpr double PI = 3.14159265358979324;
constexpr double PI2 = 2*PI;

bool is_even(int64_t n) noexcept { return n % 2 == 0; }

// range-reduces x in [0,pi] for odd n
void range_reduce_odd(double& x) noexcept
{
   if (x < 0) {
      x = -x;
   }

   if (x >= PI2) {
      x = std::fmod(x, PI2);
   }

   if (x > PI) {
      const auto p0 = 6.28125;
      const auto p1 = 0.0019353071795864769253;
      x = (p0 - x) + p1;
   }
}

// range-reduces x in [0,pi] for even n, retuns sign
double range_reduce_even(double& x) noexcept
{
   double sgn = 1.0;

   if (x < 0) {
      x = -x;
      sgn = -1;
   }

   if (x >= PI2) {
      x = std::fmod(x, PI2);
   }

   if (x > PI) {
      const auto p0 = 6.28125;
      const auto p1 = 0.0019353071795864769253;
      x = (p0 - x) + p1;
      sgn = -sgn;
   }

   return sgn;
}

// range-reduces x to be in [0,pi], returns sign
double range_reduce(int64_t n, double& x) noexcept
{
   double sgn = 1.0;

   if (is_even(n)) {
      sgn = range_reduce_even(x);
   } else {
      range_reduce_odd(x);
   }

   return sgn;
}

// returns Cl(n,0)
double Cln0(int64_t n) noexcept
{
   // zeta(n) for n = 1,...,10
   const double zeta[] = {
      1.6449340668482264, 1.2020569031595943, 1.0823232337111382,
      1.0369277551433699, 1.0173430619844491, 1.0083492773819228,
      1.0040773561979443, 1.0020083928260822, 1.0009945751278181
   };

   if (is_even(n)) {
      return 0;
   }

   return zeta[n - 2];
}

} // anonymous namespace

/**
 * @brief Clausen function \f$\operatorname{Cl}_n(\theta)\f$ for \f$n>0\f$
 * @param n degree of Clausen function
 * @param x real angle
 * @return \f$\operatorname{Cl}_n(\theta)\f$
 * @author Alexander Voigt
 */
double Cl(int64_t n, double x)
{
   if (n < 1) {
      return std::numeric_limits<double>::quiet_NaN();
   }

   if (n == 1) {
      return Cl1(x);
   }

   const auto sgn = range_reduce(n, x);

   if (x == 0 && is_even(n)) {
      return 0;
   }

   const std::complex<double> li = sgn*Li(n, std::polar(1.0, x));

   if (is_even(n)) {
      return std::imag(li);
   }

   return std::real(li);
}

} // namespace polylogarithm
