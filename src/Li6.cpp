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
   const double PI = 3.1415926535897932384626433832795;
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
   const double PI    = 3.1415926535897932384626433832795;
   const double PI2   = PI*PI;
   const double PI4   = PI2*PI2;
   const double PI6   = PI2*PI4;
   const double zeta6 = 1.0173430619844491397145179297909;
   static const int N = 40;

   const double bf[N] = {
      1., -31./64.,
      0.15241340877914951989026063100137    , -0.034365555877057613168724279835391  ,
      0.0057174797239368998628257887517147  , -0.00068180453746570644718792866941015,
      0.000049960361948734493114517016962133, -4.9166051196039047720253082216462e-7 ,
     -3.063297516130216378535303043104e-7   ,  1.4414599270849095361253742144801e-8 ,
      3.727243823092410657685992456646e-9   , -3.7300867345487607207797722915926e-10,
     -5.1246526816085832434008764611572e-11 ,  9.0541930956636682886879744762089e-12,
      6.7381882615512517077610754493801e-13 , -2.1215831150303135318527968929661e-13,
     -6.8408811719011697663920451878191e-15 ,  4.8691178462005581320638758651292e-15,
     -4.8439878499872504165088964747316e-18 , -1.1027104849107490937043030257605e-16,
      3.3353796916939381662488941184074e-18 ,  2.4735307488641352979154098748714e-18,
     -1.4370616434232492021688368713474e-19 , -5.5047110335098118061482643505979e-20,
      4.7467713917327224979130984066262e-21 ,  1.2158387178068105224373981741629e-21,
     -1.4107552403561850041424030907801e-22 , -2.6638831253268346596585643767712e-23,
      3.9667657428631007976790022608178e-24 ,  5.7821697358543615311236619348113e-25,
     -1.0787778063164257317299887699592e-25 , -1.2407397086756909899014773697759e-26,
      2.8704117917893601704260952468409e-27 ,  2.6235553563029330616574752038361e-28,
     -7.5229485465754127261588121413434e-29 , -5.4401788379624696182029193072231e-30,
      1.9502579532510166379386222338166e-30 ,  1.0978494282205187996117821359701e-31,
     -5.0149583574163009207446958541576e-32 , -2.1286737504392761053563380691739e-33
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
      const double zeta3 = 1.2020569031595942853997381615114;
      const double zeta5 = 1.0369277551433699263313654864570;
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
