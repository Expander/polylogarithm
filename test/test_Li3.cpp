#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "Li3.hpp"
#include "bench.hpp"
#include "read_data.hpp"
#include "tsil_cpp.h"
#include <cmath>
#include <utility>

#define CHECK_CLOSE(a,b,eps) CHECK((a) == doctest::Approx(b).epsilon(eps))
#define CHECK_CLOSE_COMPLEX(a,b,eps) do {                               \
      CHECK_CLOSE(std::real(a), std::real(b), (eps));                   \
      CHECK_CLOSE(std::imag(a), std::imag(b), (eps));                   \
   } while (0)
#define CHECK_SMALL(a,eps) CHECK(std::abs(a) < (eps))

template <class T> T sqr(T x) { return x*x; }
template <class T> T pow3(T x) { return x*x*x; }

std::complex<double> clog(std::complex<double> z) {
   std::complex<double> zf(z);
   // convert -0.0 to 0.0
   if (std::real(zf) == 0.0) zf.real(0.0);
   if (std::imag(zf) == 0.0) zf.imag(0.0);
   return std::log(zf);
}

std::complex<double> to_c64(std::complex<long double> z)
{
   return std::complex<double>(std::real(z), std::imag(z));
}

std::complex<long double> tsil_Li3(std::complex<long double> z) {
   return TSIL_Trilog_(z);
}

const auto Relation_1 = [](std::complex<double> z) {
   using polylogarithm::Li3;
   return Li3(z) + Li3(-z) - Li3(z*z)/4.;
};

const auto Relation_2 = [](std::complex<double> z) -> std::complex<double> {
   using polylogarithm::Li3;
   using std::log;

   if (std::abs(z) == 0)
      return {0.,0.};

   // Relation does not seem to hold for Li[3,z], not even for
   // Mathematica's PolyLog[3,z], when 0 < Re[z] < 1.
   if (std::real(z) > 0. && std::real(z) < 1.)
      return {0.,0.};

   return Li3(z) - Li3(1./z) - (-pow3(clog(-z))/6. - M_PI*M_PI/6.*clog(-z));
};

const auto Relation_3 = [](std::complex<double> z) -> std::complex<double> {
   using polylogarithm::Li3;
   using std::log;
   const double zeta3 = 1.202056903159594;

   if (std::abs(std::real(1. - z)) < 1e-10)
      return {0.,0.};

   // Relation does not seem to hold for Li[3,z], not even for
   // Mathematica's PolyLog[3,z], when Re[z] < 0 and Im[z] = 0.
   if (std::real(z) <= 0. && std::imag(z) == 0.)
      return {0.,0.};

   return Li3(z) + Li3(1.-z) + Li3(1.-1./z)
      - (zeta3 + pow3(log(z))/6. + M_PI*M_PI/6.*clog(z) - 0.5*sqr(clog(z))*clog(1.-z));
};

TEST_CASE("test_special_values")
{
   using polylogarithm::Li3;

   const double ln2 = std::log(2.);
   const double pi2 = sqr(M_PI);
   const double zeta3 = 1.2020569031595942853997381615114;
   const double phi = 0.5*(std::sqrt(5.) + 1); // golden ratio
   const double eps = 1e-15;

   const std::complex<double> zero(0.0, 0.0);
   CHECK_CLOSE_COMPLEX(Li3(zero), 0.0, eps);

   const std::complex<double> one(1.0, 0.0);
   CHECK_CLOSE_COMPLEX(Li3(one), zeta3, eps);

   const std::complex<double> mone(-1.0, 0.0);
   CHECK_CLOSE_COMPLEX(Li3(mone), -3./4.*zeta3, eps);

   const std::complex<double> half(0.5, 0.0);
   CHECK_CLOSE_COMPLEX(Li3(half), pow3(ln2)/6. - pi2/12.*ln2 + 7./8.*zeta3, eps);

   const std::complex<double> gr(1/sqr(phi), 0.0);
   CHECK_CLOSE_COMPLEX(Li3(gr), 4./5.*zeta3 + 2./3.*pow3(std::log(phi)) - 2./15.*pi2*std::log(phi), eps);
}

TEST_CASE("test_fixed_values")
{
   const auto eps64  = 1e-14;  // @todo increase precision
   const auto eps128 = 1e-15L; // @todo increase precision
   const std::string filename(std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "Li3.txt");
   const auto fixed_values = polylogarithm::test::read_from_file<long double>(filename);

   for (auto v: fixed_values) {
      const auto z128 = v.first;
      const auto z64 = to_c64(z128);
      const auto li128_expected = v.second;
      const auto li64_expected = to_c64(li128_expected);

      const auto li64_cmpl  = polylogarithm::Li3(z64);
      const auto li128_tsil = tsil_Li3(z128);

      INFO("z(128)        = " << z128);
      INFO("Li3(64)  cmpl = " << li64_expected << " (expected)");
      INFO("Li3(128) cmpl = " << li128_expected << " (expected)");
      INFO("Li3(128) cmpl = " << li128_tsil << " (TSIL)");
      INFO("Li3(64)  cmpl = " << li64_cmpl << " (polylogarithm)");

      CHECK_CLOSE_COMPLEX(li64_cmpl , li64_expected , eps64);
      CHECK_CLOSE_COMPLEX(li128_tsil, li128_expected, eps128);

      CHECK_SMALL(Relation_1(z64), eps64);
      CHECK_SMALL(Relation_2(z64), eps64);
      CHECK_SMALL(Relation_3(z64), eps64);
   }
}

TEST_CASE("test_complex_random_values")
{
   using namespace polylogarithm::bench;

   const auto values = generate_random_complexes<double>(10000, -10, 10);

   for (auto v: values) {
      const std::complex<double> li3 = polylogarithm::Li3(v);
      const std::complex<double> li3_tsil = to_c64(tsil_Li3(v));

      CHECK_CLOSE_COMPLEX(li3, li3_tsil, 1e-14);
   }
}
