#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "harmonic.hpp"

TEST_CASE("test_fixed_values")
{
   using polylogarithm::harmonic;

   CHECK(harmonic(1) == 1.0);
   CHECK(harmonic(2) == 1.5);
   CHECK(harmonic(19) == 3.5477396571436819);
   CHECK(harmonic(20) == 3.5977396571436819);
   CHECK(harmonic(100) == 5.1873775176396203);
}
