#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "zeta.hpp"
#include <cmath>

TEST_CASE("test_fixed_values")
{
   using polylogarithm::zeta;
   const double PI = 3.1415926535897932;

   CHECK(std::isinf(zeta(-263)));
   CHECK(zeta(-262) == 0.0);
   CHECK(std::isinf(zeta(-261)));
   CHECK(zeta(-259) == 8.7601563446229215e306);
   CHECK(zeta(-257) == -5.1754977470366798e303);
   CHECK(zeta(-255) == 3.1055517596048927e300);
   CHECK(zeta(-4) == 0.0);
   CHECK(zeta(-3) == 1.0/120.0);
   CHECK(zeta(-2) == 0.0);
   CHECK(zeta(-1) == -1.0/12.0);
   CHECK(zeta(0) == -0.5);
   CHECK(std::isinf(zeta(1)));
   CHECK(zeta(2) == PI*PI/6.0);
   CHECK(zeta(32) == 1.0000000002328312);
   CHECK(zeta(33) == 1.0000000001164155);
   CHECK(zeta(34) == 1.0000000000582077);
   CHECK(zeta(35) == 1.0000000000291039);
}
