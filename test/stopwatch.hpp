// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the GNU Lesser General Public
// License (GNU LGPL) version 3.
// ====================================================================

#pragma once

#include <chrono>

namespace polylogarithm {

class Stopwatch {
public:
   void start() {
      start_point = std::chrono::high_resolution_clock::now();
   }

   void stop() {
      stop_point = std::chrono::high_resolution_clock::now();
   }

   double get_time_in_seconds() const {
      typedef std::chrono::duration<int,std::micro> microseconds_t;
      microseconds_t duration(std::chrono::duration_cast<microseconds_t>(
                                 stop_point - start_point));
      return duration.count() * 0.000001;

   }

private:
   std::chrono::high_resolution_clock::time_point start_point, stop_point;
};

} // namespace polylogarithm
