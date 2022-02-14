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

   /// expansion order
   const int64_t N = 50;

   /// Bernoulli numbers B0, ..., B49
   /// Table[BernoulliB[n], {n,0,49}]
   const double bernoulli[N] = {
      1, -0.5               , 1./6.                   , 0,
      -3.333333333333333e-02, 0, 2.380952380952381e-02, 0,
      -3.333333333333333e-02, 0, 7.575757575757576e-02, 0,
      -2.531135531135531e-01, 0, 1.166666666666667e+00, 0,
      -7.092156862745098e+00, 0, 5.497117794486215e+01, 0,
      -5.291242424242424e+02, 0, 6.192123188405797e+03, 0,
      -8.658025311355312e+04, 0, 1.425517166666667e+06, 0,
      -2.729823106781609e+07, 0, 6.015808739006424e+08, 0,
      -1.511631576709215e+10, 0, 4.296146430611667e+11, 0,
      -1.371165520508833e+13, 0, 4.883323189735932e+14, 0,
      -1.929657934194006e+16, 0, 8.416930475736827e+17, 0,
      -4.033807185405945e+19, 0, 2.115074863808199e+21, 0,
      -1.208662652229652e+23, 0
   };

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

   /*

   /// Bernoulli polynomial
   std::complex<double> bernoulli_p(int64_t m, const std::complex<double>& z)
   {
      std::complex<double> result(0.0, 0.0);

      // pre-compute powers
      std::vector<std::complex<double>> powers(m + 1);
      for (int64_t k = 0; k <= m; k++) {
         powers[k] = std::pow(z + static_cast<double>(k), m);
      }

      for (int64_t n = 0; n <= m; n++) {
         std::complex<double> sum(0.0, 0.0);
         for (int64_t k = 0; k <= n; k++) {
            const double sgn = is_even(k) ? 1. : -1.;
            sum += sgn*binomial(n,k)*powers[k];
         }
         result += sum/(n + 1.);
      }

      return result;
   }

   */

   /// calculates X(p,n) for all possible n < N, p >= 0
   std::array<double,N> Xn(int64_t p)
   {
      using TArray = std::array<double,N>;
      std::vector<TArray> xn(p+1);

      // calculate X(0,n) for p = 0
      {
         TArray ar;
         for (int64_t ni = 0; ni < N; ni++) {
            ar[ni] = bernoulli[ni];
         }
         xn[0] = ar;
      }

      // pre-computing binomial coefficients
      std::array<TArray,N> binomi{};
      for (int64_t ni = 0; ni < N; ni++) {
         for (int64_t k = 0; k <= ni; k++) {
            binomi[ni][k] = binomial(ni,k);
         }
      }

      for (int64_t pi = 1; pi <= p; pi++) {
         // calculate X(pi,n) for all n < N
         TArray ar;
         for (int64_t ni = 0; ni < N; ni++) {
            double sum = 0.;
            for (int64_t k = 0; k <= ni; k++) {
               sum += binomi[ni][k]*bernoulli[ni-k]/(k+1)*xn[pi-1][k];
            }
            ar[ni] = sum;
         }
         xn[pi] = ar;
      }

      return xn[p];
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
      const std::complex<double> mu = clog(z);
      std::complex<double> sum(0.,0.), sum_old(0.,0.);
      int64_t k = 0;

      do {
         if (k == n-1) {
            k++;
            continue;
         }
         sum_old = sum;
         sum += zeta(n-k)*inv_fac(k)*std::pow(mu,k);
         k++;
      } while (!is_close(sum, sum_old, eps_d) &&
               k < std::numeric_limits<int64_t>::max() - 2);

      return std::pow(mu, n - 1)*inv_fac(n - 1)*(harmonic(n - 1) - clog(-mu)) + sum;
   }

   /// returns remainder from inversion formula
   std::complex<double> Li_rest(int64_t n, const std::complex<double>& z) noexcept
   {
      const double PI = 3.141592653589793;
      const std::complex<double> IPI2(0.,2*PI);
      const std::complex<double> lnz = clog(-z);

      // return -std::pow(IPI2, n)*inv_fac(n)*bernoulli_p(n, 0.5 + lnz/IPI2);

      std::complex<double> sum(0.0, 0.0);

      for (int64_t k = n/2; k != 0; k--) {
         const double ifac = inv_fac(n - 2*k);
         if (ifac == 0) { return 2.0*sum; }
         sum += neg_eta(2*k)*ifac*std::pow(lnz, n - 2*k);
      }

      return 2.0*sum - std::pow(lnz, n)*inv_fac(n);
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
   }

   // transformation z to y in the unit circle
   std::complex<double> y(z), r(0.0, 0.0);
   double sgn = 1;

   if (std::norm(z) > 1.0) {
      y = 1.0/z;
      r = Li_rest(n, z);
      sgn = is_even(n) ? -1. : 1.;
   }

   if (n >= 12) {
      // const auto li = Li_naive_sum(n, z);
      // const auto liy = Li_naive_sum(n, y);
      // std::cout << "n = " << n << ", z = " << z << ", y = " << y << ": Li(n,z) = " << li << ", sgn = " << sgn << ", r = " << r << ", Li(n,y) = " << liy << '\n';
      return sgn*Li_naive_sum(n, y) + r;
   }

   if (is_close(y, 1., 2e-2)) {
      return sgn*Li_expand_around_unity(n, y) + r;
   }

   const std::complex<double> u = -clog(1. - y);
   const std::array<double,N> xn = Xn(n-2);
   std::complex<double> sum(0.0, 0.0);

   for (int64_t k = N - 1; k >= 0; k--) {
      sum = u*(sum + xn[k]*inv_fac(k + 1));
   }

   return sgn*sum + r;
}

} // namespace polylogarithm
