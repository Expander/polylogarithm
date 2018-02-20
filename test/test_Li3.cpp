#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "dilog.hpp"
#include <cmath>

#define CHECK_CLOSE(a,b,eps) CHECK((a) == doctest::Approx(b).epsilon(eps))
#define CHECK_SMALL(a,eps) CHECK(std::abs(a) < (eps))
#define CHECK_CLOSE_COMPLEX(a,b,eps) do {                               \
      CHECK(std::real(a) == doctest::Approx(std::real(b)).epsilon(eps)); \
      CHECK(std::imag(a) == doctest::Approx(std::imag(b)).epsilon(eps)); \
   } while (0);

template <class T> T sqr(T x) { return x*x; }
template <class T> T pow3(T x) { return x*x*x; }

TEST_CASE("test_special_values")
{
   using namespace dilogarithm;

   const double ln2 = std::log(2.);
   const double pi2 = sqr(M_PI);
   const double zeta3 = 1.2020569031595942853997381615114;
   const double phi = 0.5*(std::sqrt(5.) + 1); // golden ratio

   CHECK_CLOSE_COMPLEX(Li3(0), 0, 1e-15);
   CHECK_CLOSE_COMPLEX(Li3(1), zeta3, 1e-15);
   CHECK_CLOSE_COMPLEX(Li3(-1), -3./4.*zeta3, 1e-15);
   CHECK_CLOSE_COMPLEX(Li3(0.5), pow3(ln2)/6. - pi2/12.*ln2 + 7./8.*zeta3, 1e-15);
   CHECK_CLOSE_COMPLEX(Li3(1/sqr(phi)), 4./5.*zeta3 + 2./3.*pow3(std::log(phi)) - 2./15.*pi2*std::log(phi), 1e-15);
}
