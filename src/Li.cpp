// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the GNU Lesser General Public
// License (GNU LGPL) version 3.
// ====================================================================

#include "Li.hpp"
#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace polylogarithm {

namespace {
   const double epsilon = std::pow(10., -std::floor(std::numeric_limits<double>::digits10));
   const double inf = std::numeric_limits<double>::infinity();

   /// expansion order
   static const long N = 50;

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

   /// 1/n! for n = 1, ..., 50
   /// Table[1/Factorial[n], {n,1,50}]
   const double fac_inv[N] = {
      1                    , 0.5                  , 1./6.                ,
      4.166666666666666e-02, 8.333333333333333e-03, 1.388888888888889e-03,
      1.984126984126984e-04, 2.480158730158730e-05, 2.755731922398589e-06,
      2.755731922398589e-07, 2.505210838544172e-08, 2.087675698786810e-09,
      1.605904383682161e-10, 1.147074559772973e-11, 7.647163731819816e-13,
      4.779477332387385e-14, 2.811457254345521e-15, 1.561920696858623e-16,
      8.220635246624330e-18, 4.110317623312165e-19, 1.957294106339126e-20,
      8.896791392450574e-22, 3.868170170630684e-23, 1.611737571096118e-24,
      6.446950284384474e-26, 2.479596263224797e-27, 9.183689863795546e-29,
      3.279889237069838e-30, 1.130996288644772e-31, 3.769987628815905e-33,
      1.216125041553518e-34, 3.800390754854744e-36, 1.151633562077195e-37,
      3.387157535521162e-39, 9.677592958631890e-41, 2.688220266286636e-42,
      7.265460179153071e-44, 1.911963205040282e-45, 4.902469756513544e-47,
      1.225617439128386e-48, 2.989310827142405e-50, 7.117406731291439e-52,
      1.655210867742195e-53, 3.761842881232262e-55, 8.359650847182804e-57,
      1.817315401561479e-58, 3.866628513960594e-60, 8.055476070751238e-62,
      1.643974708316579e-63, 3.287949416633158e-65
   };

   bool is_close(double a, double b, double eps = epsilon)
   {
      return std::abs(a - b) < eps;
   }

   bool is_close(const std::complex<double>& a, const std::complex<double>& b,
                 double eps = epsilon)
   {
      return is_close(std::real(a), std::real(b), eps) &&
             is_close(std::imag(a), std::imag(b), eps);
   }

   bool is_even(long n) { return n % 2 == 0; }

   /// complex logarithm, converts -0.0 to 0.0
   std::complex<double> clog(std::complex<double> z) noexcept
   {
      if (std::real(z) == 0.0) z.real(0.0);
      if (std::imag(z) == 0.0) z.imag(0.0);
      return std::log(z);
   }

   /// Binomial coefficients
   /// https://www.geeksforgeeks.org/space-and-time-efficient-binomial-coefficient/
   double binomial(long n, long k)
   {
      double result = 1.;

      // (n, k) = (n, n-k)
      if (k > n - k)
         k = n - k;

      for (long i = 0; i < k; i++) {
         result *= (n - i);
         result /= (i + 1);
      }

      return result;
   }

   /// Bernoulli polynomial
   std::complex<double> bernoulli_p(long m, const std::complex<double>& z)
   {
      std::complex<double> result;

      for (long n = 0; n <= m; n++) {
         std::complex<double> sum;
         for (long k = 0; k <= n; k++) {
            const double sgn = is_even(k) ? 1. : -1.;
            sum += sgn*binomial(n,k)*std::pow(z + static_cast<double>(k), m);
         }
         result += sum/(n + 1.);
      }

      return result;
   }

   /// factorial
   double fac(double n)
   {
      double result = 1.;
      for (long i = 1; i <= n; ++i)
         result *= i;
      return result;
   }

   /// Riemann zeta function for integer s (positive or negative)
   double zeta(long s)
   {
#if __cpp_lib_math_special_functions >= 201603
      return std::riemann_zeta(s);
#else
      if (s == 1)
         return inf;

      double sum = 0., sum_old = 0.;
      long n = 0;

      do {
         sum_old = sum;

         double sub_sum = 0.;

         for (long k = 0; k <= n; k++) {
            const long sgn = is_even(k) ? 1 : -1;
            sub_sum += binomial(n,k)*sgn*std::pow(k+1,-s);
         }

         sum += sub_sum*std::pow(2.,-(n+1));
         n++;
      } while (!is_close(sum_old, sum) &&
               n < std::numeric_limits<long>::max() - 2);

      return sum/(1. - std::pow(2.,1-s));
#endif
   }

   /// Dirichlet eta function
   double eta(long n)
   {
      return (1. - std::pow(2.,1-n))*zeta(n);
   }

   /// calculates X(p,n) for all possible n < N, p >= 0
   std::vector<std::array<double,N>> Xn(long p)
   {
      using TArray = std::array<double,N>;
      std::vector<TArray> xn(p+1);

      // calculate X(0,n) for p = 0
      {
         TArray ar;
         for (long ni = 0; ni < N; ni++) {
            ar[ni] = bernoulli[ni];
         }
         xn[0] = ar;
      }

      for (long pi = 1; pi <= p; pi++) {
         // calculate X(pi,n) for all n < N
         TArray ar;
         for (long ni = 0; ni < N; ni++) {
            double sum = 0.;
            for (long k = 0; k <= ni; k++) {
               sum += binomial(ni,k)*bernoulli[ni-k]/(k+1)*xn[pi-1][k];
            }
            ar[ni] = sum;
         }
         xn[pi] = ar;
      }

      return xn;
   }

   /// series expansion of Li_n(z) for n <= 0
   std::complex<double> Li_negative(long n, const std::complex<double>& z)
   {
      if (is_close(z, {1.,0.}))
         return {inf, inf};

      std::complex<double> result;

      for (long k = 0; k <= -n; k++) {
         double sum = 0.;
         for (long j = 0; j <= k; j++) {
            const long sgn = is_even(j) ? -1 : 1;
            sum += sgn*binomial(k,j)*std::pow(j+1,-n);
         }

         result += std::pow(-z/(1. - z),k+1)*sum;
      }

      if (is_close(std::imag(z), 0.))
         result.imag(0.);

      return result;
   }

   /// Series expansion of Li_n(z) in terms of powers of z.
   /// Fast convergence for large n >= 12.
   std::complex<double> Li_naive_sum(long n, const std::complex<double>& z)
   {
      std::complex<double> sum, sum_old;
      std::complex<double> pz(1.,0.);
      long k = 0;

      do {
         k++;
         pz *= z;
         sum_old = sum;
         sum += pz/std::pow(k,n);
      } while (!is_close(sum, sum_old) &&
               k < std::numeric_limits<long>::max() - 2);

      return sum;
   }

   /// Harmonic number n
   double H(long n)
   {
      double sum = 0.;

      for (long h = 1; h <= n; h++)
         sum += 1./h;

      return sum;
   }

   /// Series expansion of Li_n(z) around z ~ 1, n > 0
   std::complex<double> Li_expand_around_unity(long n, const std::complex<double>& z)
   {
      const std::complex<double> mu = clog(z);
      std::complex<double> sum, sum_old;
      long k = 0;

      do {
         if (k == n-1) {
            k++;
            continue;
         }
         sum_old = sum;
         sum += zeta(n-k)/fac(k)*std::pow(mu,k);
         k++;
      } while (!is_close(sum, sum_old) &&
               k < std::numeric_limits<long>::max() - 2);

      return std::pow(mu,n-1)/fac(n-1)*(H(n-1) - clog(-mu)) + sum;
   }

} // anonymous namespace

/**
 * @brief Clausen function \f$\mathrm{Cl}_n(\theta)\f$
 * @param n degree of Clausen function
 * @param x real angle
 * @return \f$\mathrm{Cl}_n(\theta)\f$
 */
double Cl(long n, double x)
{
   using std::exp;
   const double PI = 3.141592653589793;
   const std::complex<double> i(0.,1.);
   const std::complex<double> li = Li(n, exp(i*x));

   while (x >= 2*PI)
      x -= 2*PI;

   while (x < 0.)
      x += 2*PI;

   if (is_even(n))
      return std::imag(li);

   return std::real(li);
}

/**
 * @brief Complex polylogarithm \f$\mathrm{Li}_n(z)\f$
 * @param n degree of the polylogarithm
 * @param z complex argument
 * @return \f$\mathrm{Li}_n(z)\f$
 */
std::complex<double> Li(long n, const std::complex<double>& z)
{
   if (n < 0)
      return Li_negative(n,z);
   if (n == 0) {
      if (is_close(z, {1.,0.}))
         return {inf, inf};
      return z/(1. - z);
   }
   if (n == 1)
      return -clog(1. - z);

   if (is_close(z, 0.))
      return 0.;
   if (is_close(z, 1.))
      return zeta(n);
   if (is_close(z, -1.))
      return -eta(n);

   if (n >= 12)
      return Li_naive_sum(n, z);

   if (is_close(z, 1., 2e-2))
      return Li_expand_around_unity(n,z);

   std::complex<double> u, r;
   double sgn = 1.;

   if (std::abs(z) <= 1.) {
      u = -clog(1. - z);
   } else { // az > 1.
      const double PI = 3.141592653589793;
      const std::complex<double> IPI2(0.,2*PI);
      const std::complex<double> lnz = clog(-z);
      u = -clog(1. - 1./z);
      r = -std::pow(IPI2, n)/fac(n)*bernoulli_p(n, 0.5 + lnz/IPI2);
      sgn = is_even(n) ? -1. : 1.;
   }

   std::complex<double> p = 1., sum;

   const auto xn = Xn(n-2);

   for (long k = 0; k < N; k++) {
      p *= u;
      sum += xn[n-2][k]*p*fac_inv[k];
   }

   return sgn*sum + r;
}

} // namespace polylogarithm
