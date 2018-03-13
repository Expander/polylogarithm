// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the GNU Lesser General Public
// License (GNU LGPL) version 3.
// ====================================================================

#include "Li4.hpp"
#include <cmath>
#include <limits>

namespace polylogarithm {

namespace {
   const double epsilon = std::pow(10., -std::floor(std::numeric_limits<double>::digits10));

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

   bool is_close(const std::complex<double>& a, const std::complex<double>& b,
                 double eps = epsilon)
   {
      return std::abs(std::real(a) - std::real(b)) < eps &&
             std::abs(std::imag(a) - std::imag(b)) < eps;
   }
} // anonymous namespace

/**
 * @brief Clausen function \f$\mathrm{Cl}_4(\theta) = \mathrm{Im}(\mathrm{Li}_4(e^{i\theta}))\f$
 * @param x real angle
 * @return \f$\mathrm{Cl}_4(\theta)\f$
 */
double Cl4(double x)
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

   return std::imag(Li4(exp(i*x)));
}

/**
 * @brief Complex polylogarithm \f$\mathrm{Li}_4(z)\f$
 * @param z complex argument
 * @return \f$\mathrm{Li}_4(z)\f$
 */
std::complex<double> Li4(const std::complex<double>& z)
{
   const double PI    = 3.1415926535897932384626433832795;
   const double PI2   = PI*PI;
   const double PI4   = PI2*PI2;
   const double zeta4 = 1.0823232337111381915160036965412;
   static const int N = 40;

   const double bf[N] = {
      1., -7./16., 0.11651234567901234567901234567901, -0.019820601851851851851851851851852,
      0.0019279320987654320987654320987654, -0.000031057098765432098765432098765432,
      -0.000015624009114857835298392473643526, 8.4851235467732066371522153835079e-7,
      2.2909616603189711445359383547004e-7, -2.1832614218526916939615352313765e-8,
      -3.8828248791720155722806620380777e-9, 5.4462921032203321182579858808232e-10,
      6.9608052106827254078772334134121e-11, -1.3375737686445215199578072203635e-11,
      -1.2784852685266571604146246361574e-12, 3.2605628580248922428788418178217e-13,
      2.3647571168618257362309504812439e-14, -7.9231351220311617024299900711372e-15,
      -4.3452915709984187250497371626475e-16, 1.9236270062535920116126875526753e-16,
      7.8124143331959546707222938968737e-18, -4.6718038448036555203176282428722e-18,
      -1.3435344329812847856260222675894e-19, 1.1356826851347343244764698375938e-19,
      2.1152756202432586847505983414192e-21, -2.7642026334746517388281729253731e-21,
      -2.7068176608240064256109059558195e-23, 6.7372044828628572143267161265626e-23,
      1.3287265456683822975818009045013e-25, -1.6443773056367826467816763114889e-24,
      8.2836058999339341109829673400349e-27, 4.0190848495069350699709315007621e-26,
      -4.5757138444848790382359734346537e-28, -9.8364109094615127758320974982117e-28,
      1.6900339556037851067729523121903e-29, 2.4104805563059808504664904164902e-29,
      -5.4266127056714182501325034058929e-31, -5.9142429588741767864337599966928e-31,
      1.6232110901087370772711176143968e-32, 1.4527595437740275946132587316158e-32
   };

   if (is_close(z, 0.))
      return 0.;
   if (is_close(z, 1.))
      return zeta4;
   if (is_close(z, -1.))
      return -7.*PI4/720.;

   std::complex<double> u, r;
   double sgn = 1;

   if (std::abs(z) <= 1.) {
      u = -clog(1. - z);
   } else { // az > 1.
      const std::complex<double> lnz  = clog(-z);
      const std::complex<double> lnz2 = sqr(lnz);
      const std::complex<double> lnz4 = pow4(lnz);
      u = -clog(1. - 1./z);
      r = 1./360.*(-7*PI4 - 30.*PI2*lnz2 - 15.*lnz4);
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
