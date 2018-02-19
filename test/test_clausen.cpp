#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "dilog.hpp"
#include <cmath>

#define CHECK_CLOSE(a,b,eps) CHECK((a) == doctest::Approx(b).epsilon(eps))

TEST_CASE("test_special_values")
{
   using namespace dilogarithm;
   const double pi  = M_PI;
   const double catalan = 0.91596559417721901505460351493238411077414937428167;

   CHECK_CLOSE(clausen_2(pi/2.), catalan, 1e-15);
}
