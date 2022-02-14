// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the MIT License.
// ====================================================================

#include "Li.hpp"
#include "eta.hpp"
#include "harmonic.hpp"
#include "inv_fac.hpp"
#include "zeta.hpp"
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>

namespace polylogarithm {

namespace {
   const double eps_d = 10.0*std::numeric_limits<double>::epsilon();
   const double inf = std::numeric_limits<double>::infinity();

   bool is_close(double a, double b, double eps)
   {
      return std::abs(a - b) < eps;
   }

   bool is_close(const std::complex<double>& a, double b,
                 double eps)
   {
      return is_close(std::real(a), b, eps) &&
             is_close(std::imag(a), 0.0, eps);
   }

   bool is_close(const std::complex<double>& a, const std::complex<double>& b,
                 double eps)
   {
      return is_close(std::real(a), std::real(b), eps) &&
             is_close(std::imag(a), std::imag(b), eps);
   }

   bool is_even(int64_t n) { return n % 2 == 0; }

   /// complex logarithm, converts -0.0 to 0.0
   std::complex<double> clog(const std::complex<double>& z) noexcept
   {
      const double n = std::imag(z)*std::imag(z) + std::real(z)*std::real(z);
      double a = std::arg(z);

      if (std::imag(z) == 0.0 && a < 0.0) {
         a = -a;
      }

      return { 0.5*std::log(n), a };
   }

   /// Binomial coefficients
   /// https://www.geeksforgeeks.org/space-and-time-efficient-binomial-coefficient/
   double binomial(int64_t n, int64_t k)
   {
      double result = 1.;

      // (n, k) = (n, n-k)
      if (k > n - k) {
         k = n - k;
      }

      for (int64_t i = 0; i < k; i++) {
         result *= (n - i);
         result /= (i + 1);
      }

      return result;
   }

   // n > 0
   std::vector<double> powers_to(int64_t exponent, int64_t n)
   {
      std::vector<double> powers(n);
      powers[0] = 0.0;

      for (int64_t k = 1; k < n; k++) {
         powers[k] = std::pow(k, exponent);
      }

      return powers;
   }

   /// series expansion of Li_n(z) for n <= 0
   std::complex<double> Li_negative(int64_t n, const std::complex<double>& z)
   {
      if (is_close(z, {1.,0.}, eps_d)) {
         return {inf, inf};
      }

      const std::complex<double> frac = -z/(1. - z);
      const std::vector<double> powers = powers_to(-n, -n + 2);
      std::complex<double> result(0.,0.);

      for (int64_t k = -n; k >= 0; k--) {
         double sum = 0.;
         for (int64_t j = 0; j <= k; j++) {
            const int64_t sgn = is_even(j) ? -1 : 1;
            sum += sgn*binomial(k,j)*powers[j+1];
         }

         result = frac*(result + sum);
      }

      if (is_close(std::imag(z), 0., eps_d)) {
         result.imag(0.);
      }

      return result;
   }

   /// Series expansion of Li_n(z) in terms of powers of z.
   /// Fast convergence for large n >= 12.
   std::complex<double> Li_naive_sum(int64_t n, const std::complex<double>& z)
   {
      std::complex<double> sum(0.,0.), sum_old(0.,0.);
      std::complex<double> pz(1.,0.);
      int64_t k = 0;

      do {
         k++;
         pz *= z;
         sum_old = sum;
         sum += pz/std::pow(k,n);
      } while (!is_close(sum, sum_old, eps_d) &&
               k < std::numeric_limits<int64_t>::max() - 2);

      return sum;
   }

   /// Series expansion of Li_n(z) around z ~ 1, n > 0
   std::complex<double> Li_expand_around_unity(int64_t n, const std::complex<double>& z)
   {
      const std::complex<double> lnz = clog(z);
      const std::complex<double> lnz2 = lnz*lnz;
      std::complex<double> sum(zeta(n), 0.0), p(1.0, 0.0);

      for (int64_t j = 1; j < n - 1; ++j) {
         p *= lnz/static_cast<double>(j);
         sum += zeta(n - j)*p;
      }

      p *= lnz/static_cast<double>(n - 1);
      sum += (harmonic(n - 1) - clog(-lnz))*p;

      p *= lnz/static_cast<double>(n);
      sum += zeta(0)*p;

      p *= lnz/static_cast<double>(n + 1);
      sum += zeta(-1)*p;

      for (int64_t j = (n + 3); j < std::numeric_limits<int64_t>::max() - 2; j += 2) {
         p *= lnz2/static_cast<double>((j - 1)*j);
         const auto old_sum = sum;
         sum += zeta(n - j)*p;
         if (sum == old_sum) { break; }
      }

      return sum;
   }

   /// returns remainder from inversion formula
   std::complex<double> Li_rest(int64_t n, const std::complex<double>& z) noexcept
   {
      const double PI = 3.141592653589793;
      const std::complex<double> lnz = clog(-z);
      const std::complex<double> lnz2 = lnz*lnz;

      std::complex<double> sum(0.0, 0.0);
      std::complex<double> p(1.0, 0.0);

      if (is_even(n)) {
         for (int64_t k = n/2; k != 0; k--) {
            const double ifac = inv_fac(n - 2*k);
            if (ifac == 0) { return 2.0*sum; }
            sum += neg_eta(2*k)*ifac*p;
            p *= lnz2;
         }
      } else {
         p = lnz;
         for (int64_t k = (n - 1)/2; k != 0; k--) {
            const double ifac = inv_fac(n - 2*k);
            if (ifac == 0) { return 2.0*sum; }
            sum += neg_eta(2*k)*ifac*p;
            p *= lnz2;
         }
      }

      return 2.0*sum - p*inv_fac(n);
   }

} // anonymous namespace

/**
 * @brief Complex polylogarithm \f$\operatorname{Li}_n(z)\f$
 * @param n degree of the polylogarithm
 * @param z complex argument
 * @return \f$\operatorname{Li}_n(z)\f$
 */
std::complex<double> Li(int64_t n, const std::complex<double>& z)
{
   if (n < 0) {
      return Li_negative(n,z);
   } else if (n == 0) {
      if (is_close(z, {1.0, 0.0}, eps_d)) {
         return {inf, inf};
      }
      return z/(1.0 - z);
   } else if (n == 1) {
      return -clog(1.0 - z);
   } else if (is_close(z, 0.0, eps_d)) {
      return {0.0, 0.0};
   } else if (is_close(z, 1.0, eps_d)) {
      return {zeta(n), 0.0};
   } else if (is_close(z, -1.0, eps_d)) {
      return {neg_eta(n), 0.0};
   } else if (std::abs(z) <= 0.75) {
      return Li_naive_sum(n, z);
   } else if (std::abs(z) >= 1.4) {
      const double sgn = is_even(n) ? -1.0 : 1.0;
      return sgn*Li_naive_sum(n, 1.0/z) + Li_rest(n, z);
   }
   return Li_expand_around_unity(n, z);
}

} // namespace polylogarithm
