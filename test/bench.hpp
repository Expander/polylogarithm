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

std::vector<double> generate_random_doubles(int n, double start, double stop)
{
   std::minstd_rand gen;
   std::uniform_real_distribution<double> dist(start, stop);

   std::vector<double> v(n);
   std::generate(begin(v), end(v),
                 [&dist,&gen](){ return dist(gen); });

   return v;
}

std::vector<std::complex<double>> generate_random_complexes(
   int n, double start, double stop)
{
   const auto reals = generate_random_doubles(n, start, stop);
   const auto imags = generate_random_doubles(n, start, stop);

   std::vector<std::complex<double>> v(n);

   for (int i = 0; i < n; i++)
      v[i] = std::complex<double>(reals[i], imags[i]);

   return v;
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
