// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the MIT License.
// ====================================================================

#include "Li.hpp"
#include "eta.hpp"
#include "factorial.hpp"
#include "harmonic.hpp"
#include "zeta.hpp"
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>

namespace polylogarithm {

namespace {
   constexpr double inf = std::numeric_limits<double>::infinity();
   constexpr double nan = std::numeric_limits<double>::quiet_NaN();
   constexpr double PI = 3.1415926535897932;

   constexpr bool is_even(int64_t n) noexcept { return n % 2 == 0; }

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

   /// Series expansion of Li_n(z) in terms of powers of z.
   /// Fast convergence for large n >= 12.
   std::complex<double> Li_series(int64_t n, const std::complex<double>& z) noexcept
   {
      std::complex<double> sum(0.0, 0.0), sum_old(0.0, 0.0), p = z;
      int64_t k = 0;

      do {
         k++;
         sum_old = sum;
         sum += p;
         p *= z*std::pow(k/(1.0 + k), n);
      } while (sum != sum_old &&
               k < std::numeric_limits<int64_t>::max() - 2);

      return sum;
   }

   /// Series expansion of Li_n(z) around z ~ 1, n > 0
   std::complex<double> Li_unity_pos(int64_t n, const std::complex<double>& z) noexcept
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

   /// Series expansion of Li_n(z) around z ~ 1, n < 0
   std::complex<double> Li_unity_neg(int64_t n, const std::complex<double>& z) noexcept
   {
      const std::complex<double> lnz = clog(z);
      const std::complex<double> lnz2 = lnz*lnz;
      std::complex<double> sum = std::tgamma(1 - n)*std::pow(-lnz, n - 1);
      std::complex<double> lnzk, sum_old;
      int64_t k;

      if (is_even(n)) {
         lnzk = lnz;
         k = 1;
      } else {
         lnzk = lnz2;
         sum += zeta(n);
         k = 2;
      }

      do {
         sum_old = sum;
         sum += zeta(n - k)*inv_fac(k)*lnzk;
         lnzk *= lnz2;
         k += 2;
      } while (sum != sum_old);

      return sum;
   }

   /// returns remainder from inversion formula
   std::complex<double> Li_rest(int64_t n, const std::complex<double>& z) noexcept
   {
      if (n < 0) {
         return 0.0;
      }

      const std::complex<double> lnz = clog(-z);
      const std::complex<double> lnz2 = lnz*lnz;
      const int64_t kmax = is_even(n) ? n/2 : (n - 1)/2;
      std::complex<double> p = is_even(n) ? 1.0 : lnz;
      std::complex<double> sum(0.0, 0.0);

      for (int64_t k = kmax; k != 0; --k) {
         const double ifac = inv_fac(n - 2*k);
         if (ifac == 0) { return 2.0*sum; }
         sum += neg_eta(2*k)*ifac*p;
         p *= lnz2;
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
   if (std::isnan(std::real(z)) || std::isnan(std::imag(z))) {
      return {nan, nan};
   } else if (std::isinf(std::real(z)) || std::isinf(std::imag(z))) {
      return {-inf, 0.0};
   } else if (z == 0.0) {
      return {0.0, 0.0};
   } else if (z == 1.0) {
      if (n <= 0) {
         return {inf, inf};
      }
      return {zeta(n), 0.0};
   } else if (z == -1.0) {
      return {neg_eta(n), 0.0};
   } else if (n == -1) {
      return z/((1.0 - z)*(1.0 - z));
   } else if (n == 0) {
      return z/(1.0 - z);
   } else if (n == 1) {
      return -clog(1.0 - z);
   } else if (n < 0 && std::norm(clog(z)) < 0.512*0.512*4*PI*PI) {
      return Li_unity_neg(n, z);
   } else if (std::abs(z) <= 0.75) {
      return Li_series(n, z);
   } else if (std::abs(z) >= 1.4) {
      const double sgn = is_even(n) ? -1.0 : 1.0;
      return sgn*Li_series(n, 1.0/z) + Li_rest(n, z);
   } else if (n < 0) {
      return Li_unity_neg(n, z);
   }
   return Li_unity_pos(n, z);
}

} // namespace polylogarithm
