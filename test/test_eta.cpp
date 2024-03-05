#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "eta.hpp"
#include "test.hpp"
#include <cmath>

TEST_CASE("test_fixed_values")
{
   using polylogarithm::neg_eta;

   if (has_inf()) {
      CHECK(std::isinf(neg_eta(-221)));
      CHECK(std::isinf(neg_eta(-219)));
   }
   CHECK(neg_eta(-222) == 0.0);
   CHECK(neg_eta(-220) == 0.0);
   CHECK(neg_eta(-218) == 0.0);
   CHECK(neg_eta(-110) == 0.0);
   CHECK(neg_eta(-3) == 1.0/8);
   CHECK(neg_eta(-2) == 0.0);
   CHECK(neg_eta(-1) == -0.25);
   CHECK(neg_eta(0)  == -0.5);
   CHECK(neg_eta(1)  == -0.69314718055994531);
   CHECK(neg_eta(2)  == -0.82246703342411322);
   CHECK(neg_eta(52) == -0.9999999999999998);
   CHECK(neg_eta(53) == -0.9999999999999999);
   CHECK(neg_eta(54) == -0.9999999999999999);
   CHECK(neg_eta(55) == -1.0);
   CHECK(neg_eta(56) == -1.0);
}
