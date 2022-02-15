#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "factorial.hpp"

TEST_CASE("test_fixed_values")
{
   using polylogarithm::inv_fac;

   CHECK(inv_fac(  0) == 1.0);
   CHECK(inv_fac(  1) == 1.0);
   CHECK(inv_fac(  2) == 0.5);
   CHECK(inv_fac(177) == 2.8547896502574379e-323);
   CHECK(inv_fac(178) == 0.0);
}
