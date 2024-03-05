// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the MIT License.
// ====================================================================

#pragma once

#include <cmath>
#include <complex>

namespace polylogarithm {

template <typename T>
struct Complex {
   constexpr Complex(T re_ = T{}, T im_ = T{}) : re(re_), im(im_) {}
   operator std::complex<T>() const noexcept { return std::complex<T>(re, im); }
   T re{};
   T im{};
};

template <typename T>
constexpr T arg(const Complex<T>& z) noexcept
{
   return std::atan2(z.im, z.re);
}

template <typename T>
constexpr Complex<T> conj(const Complex<T>& z) noexcept
{
   return { z.re, -z.im };
}

template <typename T>
Complex<T> log(const Complex<T>& z) noexcept
{
   if (z.im == T(0) && z.re > T(0)) {
      return { std::log(z.re), T(0) };
   } else if (z.im == T(0)) {
      return { std::log(-z.re), 4*std::atan(T(1)) };
   }
   return { std::log(norm(z)), arg(z) };
}

template <typename T>
Complex<T> log1p(const Complex<T>& z) noexcept
{
   const Complex<T> u = T(1) + z;

   if (u.re == T(1) && u.im == T(0)) {
      return z;
   } else if (u.re <= T(0)) {
      return log(u);
   }

   return log(u)*(z/(u - T(1)));
}

template <typename T>
constexpr T norm(const Complex<T>& z) noexcept
{
   return std::hypot(z.re, z.im);
}

template <typename T>
constexpr T norm_sqr(const Complex<T>& z) noexcept
{
   return z.re*z.re + z.im*z.im;
}

template <typename T>
constexpr Complex<T> operator+(const Complex<T>& a, const Complex<T>& b) noexcept
{
   return { a.re + b.re, a.im + b.im };
}

template <typename T>
constexpr Complex<T> operator+(const Complex<T>& z, T x) noexcept
{
   return { z.re + x, z.im };
}

template <typename T>
constexpr Complex<T> operator+(T x, const Complex<T>& z) noexcept
{
   return { x + z.re, z.im };
}

template <typename T>
constexpr Complex<T> operator-(const Complex<T>& a, const Complex<T>& b) noexcept
{
   return { a.re - b.re, a.im - b.im };
}

template <typename T>
constexpr Complex<T> operator-(T x, const Complex<T>& z) noexcept
{
   return { x - z.re, -z.im };
}

template <typename T>
constexpr Complex<T> operator-(const Complex<T>& z, T x) noexcept
{
   return { z.re - x, z.im };
}

template <typename T>
constexpr Complex<T> operator-(const Complex<T>& z) noexcept
{
   return { -z.re, -z.im };
}

template <typename T>
constexpr Complex<T> operator*(const Complex<T>& a, const Complex<T>& b) noexcept
{
   return { a.re*b.re - a.im*b.im, a.re*b.im + a.im*b.re };
}

template <typename T>
constexpr Complex<T> operator*(T x, const Complex<T>& z) noexcept
{
   return { x*z.re, x*z.im };
}

template <typename T>
constexpr Complex<T> operator*(const Complex<T>& z, T x) noexcept
{
   return x*z;
}

template <typename T>
constexpr Complex<T> operator/(T x, const Complex<T>& z) noexcept
{
   return x*conj(z)/norm_sqr(z);
}

template <typename T>
constexpr Complex<T> operator/(const Complex<T>& z, T x) noexcept
{
   return { z.re/x, z.im/x };
}

template <typename T>
Complex<T> operator/(const Complex<T>& z, const Complex<T>& w) noexcept
{
   const T a = z.re, b = z.im, c = w.re, d = w.im;

   if (std::abs(c) >= std::abs(d)) {
      const T r = d/c;
      const T den = c + d*r;
      return { (a + b*r)/den, (b - a*r)/den };
   }

   const T r = c/d;
   const T den = c*r + d;

   return { (a*r + b)/den, (b*r - a)/den };
}

} // namespace polylogarithm
