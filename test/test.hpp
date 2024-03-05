#pragma once

#include <cmath>
#include <complex>


#ifndef M_PI
#define M_PI 3.1415926535897932
#endif


#define CHECK_CLOSE(a,b,eps) do {                       \
      if (std::isinf(a) && std::isinf(b))               \
         CHECK(true);                                   \
      else                                              \
         CHECK((a) == doctest::Approx(b).epsilon(eps)); \
   } while (0);


#define CHECK_CLOSE_COMPLEX(a,b,eps) do {               \
      CHECK_CLOSE(std::real(a), std::real(b), (eps));   \
      CHECK_CLOSE(std::imag(a), std::imag(b), (eps));   \
   } while (0)


#define CHECK_SMALL(a,eps) CHECK(std::abs(a) < (eps))


inline bool has_signed_zero() noexcept
{
   const auto fn = [] (double x) { return x == 0.0 ? x : x; };
   return std::signbit(fn(-0.0));
}


template <typename T, typename U>
std::complex<T> to(std::complex<U> z)
{
   return std::complex<T>(static_cast<T>(std::real(z)), static_cast<T>(std::imag(z)));
}
