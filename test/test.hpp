#pragma once

#include <cmath>
#include <complex>


#ifndef M_PI
#define M_PI 3.1415926535897932
#endif


#define CHECK_CLOSE(a,b,eps)                            \
   do {                                                 \
      if (std::isinf(a) && std::isinf(b)) {             \
         CHECK(true);                                   \
      } else {                                          \
         CHECK((a) == doctest::Approx(b).epsilon(eps)); \
      }                                                 \
   } while (0);


#define CHECK_CLOSE_COMPLEX(a,b,eps) do {               \
      CHECK_CLOSE(std::real(a), std::real(b), (eps));   \
      CHECK_CLOSE(std::imag(a), std::imag(b), (eps));   \
   } while (0)


#define CHECK_SMALL(a,eps) CHECK(std::abs(a) < (eps))


#define CHECK_CLOSE_REL(a,b,eps) do {                                   \
      const bool pred = is_close_rel((a), (b), (eps));                  \
      INFO("Comparing numbers " << std::setprecision(17) << (a) << " =?= " << (b) << " with relative precision " << (eps)); \
      CHECK(pred);                                                      \
   } while (0)


inline bool has_inf() noexcept
{
   return std::isinf(std::numeric_limits<double>::max() + 1.0);
}


inline bool has_signed_zero() noexcept
{
   const auto fn = [] (double x) { return x == 0.0 ? x : x; };
   return std::signbit(fn(-0.0));
}


inline bool is_ieee754_compliant() noexcept
{
   constexpr bool no_fast_math = (1.0 + 0.1) - 1.0 != 0.1;

   return has_inf() && has_signed_zero() && no_fast_math;
}


inline bool is_close_rel(double x, double y, double eps) noexcept
{
   const double ma = std::max(std::abs(x), std::abs(y));
   if (ma == 0.0) {
      return true;
   }
   return std::abs(x - y)/ma < std::abs(eps);
}


template <class T> T sqr(T x) {
   return x*x;
}


template <class T> T pow3(T x) {
   return x*x*x;
}


template <class T> T pow4(T x) {
   return sqr(sqr(x));
}


template <typename T, typename U>
std::complex<T> to(std::complex<U> z)
{
   return std::complex<T>(static_cast<T>(std::real(z)), static_cast<T>(std::imag(z)));
}
