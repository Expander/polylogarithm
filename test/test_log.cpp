#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "complex.hpp"
#include "test.hpp"
#include <complex>
#include <limits>

// tests log(z) and log1p(1+z) for small complex z
TEST_CASE("test_log_and_log1p")
{
   if (!is_ieee754_compliant()) {
      return;
   }

   using namespace polylogarithm;

   const double eps = std::numeric_limits<double>::epsilon();
   // z = 9007155738389423 / 2^53 - (5069303045819176 / 2^60) I
   const std::complex<double> z(0.999995168714454, -0.004396919500211628);
   const std::complex<double> cpp_log(std::log(z));

   const Complex<double> pl_z(std::real(z), std::imag(z));
   const Complex<double> pl_log = log(pl_z);
   const Complex<double> pl_log1p = log1p(pl_z);

   const std::complex<double> expected_log(4.8351532915892516848328240700759518475389916e-6, 0.00439691240793594643979611108925073970426108463713022280825);
   const std::complex<double> expected_log1p(0.693147181532726411790714997046464927858457, -0.002198461518912911475356949677977809228317);

   INFO("z             = " << std::setprecision(17) << z);
   INFO("std::log(z)   = " << cpp_log);
   INFO("log(z)        = " << std::complex<double>(pl_log));
   INFO("exp. log(z)   = " << expected_log);
   INFO("log1p(z)      = " << std::complex<double>(pl_log1p));
   INFO("exp. log1p(z) = " << expected_log1p);

   CHECK_CLOSE_REL(std::real(cpp_log), std::real(expected_log), eps);
   CHECK_CLOSE_REL(pl_log.re         , std::real(expected_log), eps);
   CHECK_CLOSE_REL(pl_log1p.re       , std::real(expected_log1p), eps);
}
