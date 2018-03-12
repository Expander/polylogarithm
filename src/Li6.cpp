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
   template <typename T>
   T sqr(T x) noexcept { return x*x; }

   template <typename T>
   T pow4(T x) noexcept { return x*x*x*x; }

   // converts -0.0 to 0.0
   std::complex<double> clog(std::complex<double> z) noexcept {
      if (std::real(z) == 0.0) z.real(0.0);
      if (std::imag(z) == 0.0) z.imag(0.0);
      return std::log(z);
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

   while (x > 2*PI)
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
      0.15241340877914951989026063100137, -0.034365555877057613168724279835391,
      0.0057174797239368998628257887517147, -0.00068180453746570644718792866941015,
      0.000049960361948734493114517016962133, -4.9166051196039047720253082216462e-7,
      -3.063297516130216378535303043104e-7, 1.4414599270849095361253742144801e-8,
      3.727243823092410657685992456646e-9, -3.7300867345487607207797722915926e-10,
      -5.1246526816085832434008764611572e-11, 9.0541930956636682886879744762089e-12,
      6.7381882615512517077610754493801e-13, -2.1215831150303135318527968929661e-13,
      -6.8408811719011697663920451878191e-15, 4.8691178462005581320638758651292e-15,
      -4.8439878499872504165088964747316e-18, -1.1027104849107490937043030257605e-16,
      3.3353796916939381662488941184074e-18, 2.4735307488641352979154098748714e-18,
      -1.4370616434232492021688368713474e-19, -5.5047110335098118061482643505979e-20,
      4.7467713917327224979130984066262e-21, 1.2158387178068105224373981741629e-21,
      -1.4107552403561850041424030907801e-22, -2.6638831253268346596585643767712e-23,
      3.9667657428631007976790022608178e-24, 5.7821697358543615311236619348113e-25,
      -1.0787778063164257317299887699592e-25, -1.2407397086756909899014773697759e-26,
      2.8704117917893601704260952468409e-27, 2.6235553563029330616574752038361e-28,
      -7.5229485465754127261588121413434e-29, -5.4401788379624696182029193072231e-30,
      1.9502579532510166379386222338166e-30, 1.0978494282205187996117821359701e-31,
      -5.0149583574163009207446958541576e-32, -2.1286737504392761053563380691739e-33
   };

   if (z == 0.)
      return 0.;
   if (z == 1.)
      return zeta6;
   if (z == -1.)
      return -31.*zeta6/32.;

   std::complex<double> u, r;
   double sgn = 1;

   if (std::abs(z) <= 1.) {
      u = -clog(1. - z);
   } else { // az > 1.
      const std::complex<double> lnz  = clog(-z);
      const std::complex<double> lnz2 = sqr(lnz);
      const std::complex<double> lnz4 = pow4(lnz);
      const std::complex<double> lnz6 = lnz2*lnz4;
      u = -clog(1. - 1./z);
      r = -31.*PI6/15120. - 7./720.*PI4*lnz2 - 1./144.*PI2*lnz4 - 1./720.*lnz6;
      sgn = -1;
   }

   std::complex<double> p = 1., sum;

   for (const double b: bf) {
      p *= u;
      sum += b*p;
   }

   return sgn*sum + r;
}

} // namespace polylogarithm
