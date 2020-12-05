// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the MIT License.
// ====================================================================

#pragma once

#include <chrono>

namespace polylogarithm {

class Stopwatch {
private:
   using microseconds_t = std::chrono::duration<int,std::micro>;

   microseconds_t::rep get_ticks() const {
      microseconds_t duration(std::chrono::duration_cast<microseconds_t>(
                                 stop_point - start_point));
      return duration.count();
   }

public:
   void start() {
      start_point = std::chrono::high_resolution_clock::now();
   }

   void stop() {
      stop_point = std::chrono::high_resolution_clock::now();
   }

   double get_time_in_seconds() const {
      return get_ticks() * 0.000001;
   }

   double get_time_in_milliseconds() const {
      return get_ticks() * 0.001;
   }

private:
   std::chrono::high_resolution_clock::time_point start_point, stop_point;
};

} // namespace polylogarithm
