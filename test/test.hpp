#pragma once

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
