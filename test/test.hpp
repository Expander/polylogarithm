#pragma once

inline bool has_signed_zero() noexcept
{
   const auto fn = [] (double x) { return x == 0.0 ? x : x; };
   return std::signbit(fn(-0.0));
}
