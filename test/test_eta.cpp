#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "eta.hpp"

TEST_CASE("test_fixed_values")
{
   using polylogarithm::neg_eta;

   CHECK(neg_eta(1)  == -0.69314718055994531);
   CHECK(neg_eta(2)  == -0.82246703342411322);
   CHECK(neg_eta(52) == -0.9999999999999998);
   CHECK(neg_eta(53) == -0.9999999999999999);
   CHECK(neg_eta(54) == -0.9999999999999999);
   CHECK(neg_eta(55) == -1.0);
   CHECK(neg_eta(56) == -1.0);
}
