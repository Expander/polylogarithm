// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the GNU Lesser General Public
// License (GNU LGPL) version 3.
// ====================================================================

#include "Li6.hpp"
#include <cmath>
#include <limits>

namespace polylogarithm {

namespace {
   const double epsilon = std::pow(10., -std::floor(std::numeric_limits<double>::digits10));

   template <typename T> T pow2(T x) noexcept { return x*x; }
   template <typename T> T pow4(T x) noexcept { return x*x*x*x; }

   // converts -0.0 to 0.0
   std::complex<double> clog(std::complex<double> z) noexcept {
      if (std::real(z) == 0.0) z.real(0.0);
      if (std::imag(z) == 0.0) z.imag(0.0);
      return std::log(z);
   }

   bool is_close(const std::complex<double>& a, const std::complex<double>& b,
                 double eps = epsilon)
   {
      return std::abs(std::real(a) - std::real(b)) < eps &&
             std::abs(std::imag(a) - std::imag(b)) < eps;
   }
} // anonymous namespace

/**
 * @brief Clausen function \f$\mathrm{Cl}_6(\theta) = \mathrm{Im}(\mathrm{Li}_6(e^{i\theta}))\f$
 * @param x real angle
 * @return \f$\mathrm{Cl}_6(\theta)\f$
 */
double Cl6(double x)
{
   using std::exp;
   const double PI = 3.141592653589793;
   const std::complex<double> i(0.,1.);

   while (x >= 2*PI)
      x -= 2*PI;

   while (x < 0.)
      x += 2*PI;

   if (std::abs(x) < std::numeric_limits<double>::epsilon() ||
       std::abs(x - PI) < std::numeric_limits<double>::epsilon() ||
       std::abs(x - 2*PI) < std::numeric_limits<double>::epsilon())
      return 0.;

   return std::imag(Li6(exp(i*x)));
}

/**
 * @brief Complex polylogarithm \f$\mathrm{Li}_6(z)\f$
 * @param z complex argument
 * @return \f$\mathrm{Li}_6(z)\f$
 */
std::complex<double> Li6(const std::complex<double>& z)
{
   const double PI    = 3.141592653589793;
   const double PI2   = PI*PI;
   const double PI4   = PI2*PI2;
   const double PI6   = PI2*PI4;
   const double zeta6 = 1.017343061984449;
   static const int N = 40;

   const double bf[N] = {
      1., -31./64.,
      1.524134087791495e-01, -3.436555587705761e-02,
      5.717479723936899e-03, -6.818045374657064e-04,
      4.996036194873449e-05, -4.916605119603904e-07,
     -3.063297516130216e-07,  1.441459927084909e-08,
      3.727243823092410e-09, -3.730086734548760e-10,
     -5.124652681608583e-11,  9.054193095663668e-12,
      6.738188261551251e-13, -2.121583115030313e-13,
     -6.840881171901169e-15,  4.869117846200558e-15,
     -4.843987849987250e-18, -1.102710484910749e-16,
      3.335379691693938e-18,  2.473530748864135e-18,
     -1.437061643423249e-19, -5.504711033509811e-20,
      4.746771391732722e-21,  1.215838717806810e-21,
     -1.410755240356185e-22, -2.663883125326834e-23,
      3.966765742863100e-24,  5.782169735854361e-25,
     -1.078777806316425e-25, -1.240739708675690e-26,
      2.870411791789360e-27,  2.623555356302933e-28,
     -7.522948546575412e-29, -5.440178837962469e-30,
      1.950257953251016e-30,  1.097849428220518e-31,
     -5.014958357416300e-32, -2.128673750439276e-33
   };

   if (is_close(z, 0.))
      return 0.;
   if (is_close(z, 1.))
      return zeta6;
   if (is_close(z, 1., 0.02)) {
      const std::complex<double> I(0.,1.);
      const std::complex<double> IPI(0.,PI);
      const std::complex<double> zm1 = z - 1.;
      const std::complex<double> lzm1 = clog(zm1);
      const double zeta3 = 1.202056903159594;
      const double zeta5 = 1.036927755143369;
      const double ceil = std::arg(zm1) > 0. ? 1. : 0.;
      const std::complex<double> cs[8] = {
         zeta5,
         (PI4 - 90*zeta5)/180.,
         (-PI4 + 30*(zeta3 + 2*zeta5))/180.,
         (15*PI2 + 11*PI4 - 540*(zeta3 + zeta5))/2160.,
         (-lzm1/120. + (411. + 180.*(-1 + 2*ceil)*IPI - 100*PI2*(3 + PI2)
                        + 6300*zeta3 + 4320*zeta5)/21600.),
         (lzm1/48. + (2700*(1 - 2*ceil)*IPI + 2550*PI2 + 548*PI4
                      - 45*(127 + 900*zeta3 + 480*zeta5))/129600.),
         (0.06919642857142858 + 5./144.*(-1 + 2*ceil)*IPI - 5./144.*lzm1
          - 7./7200.*PI2*(25 + 4*PI2) + 29./90.*zeta3 + zeta5/7.),
         (-0.0921875 + 7./144.*(1 - 2*ceil)*IPI + 7./144.*lzm1
                      + 967./34560.*PI2 + 121./33600.*PI4 - 469./1440.*zeta3 - zeta5/8.)
      };

      return zeta6 +
         zm1 * (cs[0] +
         zm1 * (cs[1] +
         zm1 * (cs[2] +
         zm1 * (cs[3] +
         zm1 * (cs[4] +
         zm1 * (cs[5] +
         zm1 * (cs[6] +
         zm1 * (cs[7]))))))));
   }
   if (is_close(z, -1.))
      return -31.*zeta6/32.;

   std::complex<double> u, r;
   double sgn = 1;

   if (std::abs(z) <= 1.) {
      u = -clog(1. - z);
   } else { // az > 1.
      const std::complex<double> lnz  = clog(-z);
      const std::complex<double> lnz2 = pow2(lnz);
      const std::complex<double> lnz4 = pow4(lnz);
      const std::complex<double> lnz6 = lnz2*lnz4;
      u = -clog(1. - 1./z);
      r = -31.*PI6/15120. - 7./720.*PI4*lnz2 - 1./144.*PI2*lnz4 - 1./720.*lnz6;
      sgn = -1;
   }

   const std::complex<double> sum =
      u * (bf[0] +
      u * (bf[1] +
      u * (bf[2] +
      u * (bf[3] +
      u * (bf[4] +
      u * (bf[5] +
      u * (bf[6] +
      u * (bf[7] +
      u * (bf[8] +
      u * (bf[9] +
      u * (bf[10] +
      u * (bf[11] +
      u * (bf[12] +
      u * (bf[13] +
      u * (bf[14] +
      u * (bf[15] +
      u * (bf[16] +
      u * (bf[17] +
      u * (bf[18] +
      u * (bf[19] +
      u * (bf[20] +
      u * (bf[21] +
      u * (bf[22] +
      u * (bf[23] +
      u * (bf[24] +
      u * (bf[25] +
      u * (bf[26] +
      u * (bf[27] +
      u * (bf[28] +
      u * (bf[29] +
      u * (bf[30] +
      u * (bf[31] +
      u * (bf[32] +
      u * (bf[33] +
      u * (bf[34] +
      u * (bf[35] +
      u * (bf[36] +
      u * (bf[37] +
      u * (bf[38] +
      u * (bf[39]))))))))))))))))))))))))))))))))))))))));

   return sgn*sum + r;
}

} // namespace polylogarithm
