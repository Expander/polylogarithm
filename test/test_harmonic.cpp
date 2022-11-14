#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "harmonic.hpp"

#define CHECK_CLOSE(a,b,eps) CHECK((a) == doctest::Approx(b).epsilon(eps))

TEST_CASE("test_fixed_values")
{
   using polylogarithm::harmonic;

   const double eps = 1e-15;

   CHECK_CLOSE(harmonic(1), 1.0, eps);
   CHECK_CLOSE(harmonic(2), 1.5, eps);
   CHECK_CLOSE(harmonic(19), 3.5477396571436819, eps);
   CHECK_CLOSE(harmonic(20), 3.5977396571436819, eps);
   CHECK_CLOSE(harmonic(100), 5.1873775176396203, eps);
}
