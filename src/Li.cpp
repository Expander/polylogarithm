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
      1, -0.5               , 0.16666666666666666     , 0,
      -0.03333333333333333  , 0, 0.023809523809523808 , 0,
      -0.03333333333333333  , 0, 0.07575757575757576  , 0,
      -0.2531135531135531   , 0, 1.1666666666666667   , 0,
      -7.092156862745098    , 0, 54.971177944862156   , 0,
      -529.1242424242424    , 0, 6192.123188405797    , 0,
      -86580.25311355312    , 0, 1.4255171666666667e6 , 0,
      -2.7298231067816094e7 , 0, 6.015808739006424e8  , 0,
      -1.5116315767092157e10, 0, 4.296146430611667e11 , 0,
      -1.3711655205088334e13, 0, 4.883323189735932e14 , 0,
      -1.9296579341940068e16, 0, 8.416930475736827e17 , 0,
      -4.0338071854059454e19, 0, 2.1150748638081993e21, 0,
      -1.2086626522296526e23, 0
   };

   /// 1/n! for n = 1, ..., 50
   /// Table[1/Factorial[n], {n,1,50}]
   const double fac_inv[N] = {
      1                     , 0.5                   , 0.16666666666666666    ,
      0.041666666666666664  , 0.008333333333333333  , 0.001388888888888889   ,
      0.0001984126984126984 , 0.0000248015873015873 , 2.7557319223985893e-6  ,
      2.755731922398589e-7  , 2.505210838544172e-8  , 2.08767569878681e-9    ,
      1.6059043836821613e-10, 1.1470745597729725e-11, 7.647163731819816e-13  ,
      4.779477332387385e-14 , 2.8114572543455206e-15, 1.5619206968586225e-16 ,
      8.22063524662433e-18  , 4.110317623312165e-19 , 1.9572941063391263e-20 ,
      8.896791392450574e-22 , 3.8681701706306835e-23, 1.6117375710961184e-24 ,
      6.446950284384474e-26 , 2.4795962632247972e-27, 9.183689863795546e-29  ,
      3.279889237069838e-30 , 1.1309962886447718e-31, 3.7699876288159054e-33 ,
      1.2161250415535181e-34, 3.800390754854744e-36 , 1.151633562077195e-37  ,
      3.387157535521162e-39 , 9.67759295863189e-41  , 2.688220266286636e-42  ,
      7.265460179153071e-44 , 1.911963205040282e-45 , 4.902469756513544e-47  ,
      1.225617439128386e-48 , 2.9893108271424046e-50, 7.117406731291439e-52  ,
      1.6552108677421951e-53, 3.7618428812322616e-55, 8.359650847182804e-57  ,
      1.817315401561479e-58 , 3.866628513960594e-60 , 8.055476070751238e-62  ,
      1.643974708316579e-63 , 3.2879494166331584e-65
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

      // C(n, k) = C(n, n-k)
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
      const double eps = std::pow(10., -std::floor(std::numeric_limits<double>::digits10));
      std::complex<double> sum, sum_old;
      long k = 0;

      do {
         k++;
         sum_old = sum;
         sum += std::pow(z,k)/std::pow(k,n);
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
   const double PI = 3.1415926535897932384626433832795;
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
      const double PI = 3.1415926535897932384626433832795;
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
