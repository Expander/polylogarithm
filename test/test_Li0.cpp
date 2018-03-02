#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "Li0.hpp"
#include <cmath>
#include <utility>

#define CHECK_CLOSE(a,b,eps) CHECK((a) == doctest::Approx(b).epsilon(eps))
#define CHECK_CLOSE_COMPLEX(a,b,eps) do {                               \
      CHECK(std::real(a) == doctest::Approx(std::real(b)).epsilon(eps)); \
      CHECK(std::imag(a) == doctest::Approx(std::imag(b)).epsilon(eps)); \
   } while (0);
#define CHECK_SMALL(a,eps) CHECK(std::abs(a) < (eps))

TEST_CASE("test_special_values")
{
   using polylogarithm::Li0;

   CHECK_SMALL(Li0(0.), 1e-15);
   CHECK(!std::isfinite(Li0(1)));
}
