// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the GNU Lesser General Public
// License (GNU LGPL) version 3.
// ====================================================================

#pragma once

#include "stopwatch.hpp"
#include <algorithm>
#include <complex>
#include <random>
#include <vector>

namespace polylogarithm {
namespace bench {

template <class T>
std::vector<T> generate_random_scalars(int n, T start, T stop)
{
   static std::minstd_rand gen;
   std::uniform_real_distribution<T> dist(start, stop);

   std::vector<T> v(n);
   std::generate(begin(v), end(v),
                 [&dist](){ return dist(gen); });

   return v;
}

template <class T>
std::vector<std::complex<T>> generate_random_complexes(
   int n, T start, T stop)
{
   const auto reals = generate_random_scalars<T>(n, start, stop);
   const auto imags = generate_random_scalars<T>(n, start, stop);

   std::vector<std::complex<T>> v(n);

   for (int i = 0; i < n; i++) {
      v[i] = std::complex<T>(reals[i], imags[i]);
   }

   return v;
}

template <class T>
inline void do_not_optimize(const T& value)
{
   asm volatile("" : : "r,m"(value) : "memory");
}

template <class F>
double time_in_seconds(F&& f)
{
   polylogarithm::Stopwatch sw;
   sw.start();
   f();
   sw.stop();
   return sw.get_time_in_seconds();
}

} // namespace bench
} // namespace polylogarithm
