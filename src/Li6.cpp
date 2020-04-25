// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the GNU Lesser General Public
// License (GNU LGPL) version 3.
// ====================================================================

#include "Li6.hpp"
#include <cfloat>
#include <cmath>
#include <limits>

namespace polylogarithm {

namespace {
   template <typename T> T sqr(T x) noexcept { return x*x; }

   template <typename T>
   std::complex<T> clog(std::complex<T> z) noexcept
   {
      // converts -0.0 to 0.0
      if (std::real(z) == T(0)) { z.real(T(0)); }
      if (std::imag(z) == T(0)) { z.imag(T(0)); }

      const T rz = std::real(z);
      const T iz = std::imag(z);
      const T nz = sqr(rz) + sqr(iz);

      return std::complex<T>(0.5*std::log(nz), std::atan2(iz, rz));
   }

   template <typename T>
   bool is_close(const std::complex<T>& a, T b, T eps)
   {
      return std::abs(std::real(a) - b) < eps && std::abs(std::imag(a)) < eps;
   }

   template <typename T>
   std::complex<T> cadd(T a, const std::complex<T>& b) noexcept
   {
      return std::complex<T>(a + std::real(b), std::imag(b));
   }

   template <typename T>
   std::complex<T> cadd(const std::complex<T>& a, const std::complex<T>& b) noexcept
   {
      return std::complex<T>(std::real(a) + std::real(b),
                             std::imag(a) + std::imag(b));
   }

   template <typename T>
   std::complex<T> cmul(const std::complex<T>& a, T b) noexcept
   {
      return std::complex<T>(std::real(a) * b, std::imag(a) * b);
   }

   template <typename T>
   std::complex<T> cmul(const std::complex<T>& a, const std::complex<T>& b) noexcept
   {
      return std::complex<T>(
         std::real(a) * std::real(b) - std::imag(a) * std::imag(b),
         std::real(a) * std::imag(b) + std::imag(a) * std::real(b));
   }
} // anonymous namespace

/**
 * @brief Clausen function \f$\mathrm{Cl}_6(\theta) = \mathrm{Im}(\mathrm{Li}_6(e^{i\theta}))\f$
 * @param x real angle
 * @return \f$\mathrm{Cl}_6(\theta)\f$
 */
double Cl6(double x)
{
   const double PI = 3.141592653589793;
   const std::complex<double> i(0.0, 1.0);

   while (x >= 2*PI) {
      x -= 2*PI;
   }

   while (x < 0.0) {
      x += 2*PI;
   }

   if (std::abs(x) < std::numeric_limits<double>::epsilon() ||
       std::abs(x - PI) < std::numeric_limits<double>::epsilon() ||
       std::abs(x - 2*PI) < std::numeric_limits<double>::epsilon()) {
      return 0.0;
   }

   return std::imag(Li6(std::exp(i*x)));
}

/**
 * @brief Clausen function \f$\mathrm{Cl}_6(\theta) = \mathrm{Im}(\mathrm{Li}_6(e^{i\theta}))\f$ with long double precision
 * @param x real angle
 * @return \f$\mathrm{Cl}_6(\theta)\f$
 */
long double Cl6(long double x)
{
   const long double PI = 3.14159265358979323846264338327950288L;
   const std::complex<long double> i(0.0L, 1.0L);

   while (x >= 2*PI) {
      x -= 2*PI;
   }

   while (x < 0.0L) {
      x += 2*PI;
   }

   if (std::abs(x) < std::numeric_limits<long double>::epsilon() ||
       std::abs(x - PI) < std::numeric_limits<long double>::epsilon() ||
       std::abs(x - 2*PI) < std::numeric_limits<long double>::epsilon()) {
      return 0.0L;
   }

   return std::imag(Li6(std::exp(i*x)));
}

/**
 * @brief Complex polylogarithm \f$\mathrm{Li}_6(z)\f$
 * @param z complex argument
 * @return \f$\mathrm{Li}_6(z)\f$
 */
std::complex<double> Li6(const std::complex<double>& z)
{
   const double eps   = 10.0*std::numeric_limits<double>::epsilon();
   const double PI    = 3.141592653589793;
   const double PI2   = PI*PI;
   const double PI4   = PI2*PI2;
   const double PI6   = PI2*PI4;
   const double zeta6 = 1.017343061984449;
   const double bf[18] = {
      1.0, -31.0/64.0,
      1.524134087791495e-01, -3.436555587705761e-02,
      5.717479723936899e-03, -6.818045374657064e-04,
      4.996036194873449e-05, -4.916605119603904e-07,
     -3.063297516130216e-07,  1.441459927084909e-08,
      3.727243823092410e-09, -3.730086734548760e-10,
     -5.124652681608583e-11,  9.054193095663668e-12,
      6.738188261551251e-13, -2.121583115030313e-13,
     -6.840881171901169e-15,  4.869117846200558e-15
   };

   if (is_close(z, 0.0, eps)) {
      return { 0.0, 0.0 };
   }
   if (is_close(z, 1.0, eps)) {
      return { zeta6, 0.0 };
   }
   if (is_close(z, -1.0, eps)) {
      return { -31.0*zeta6/32.0, 0.0 };
   }

   const auto az  = std::abs(z);
   const auto pz  = std::arg(z);
   const auto lnz = std::log(az);

   if (sqr(lnz) + sqr(pz) < 1.0) { // |log(z)| < 1
      const auto u  = std::complex<double>(lnz, pz); // clog(z)
      const auto u2 = u*u;
      const auto c0 = zeta6;
      const auto c1 = 1.036927755143370; // zeta(5)
      const auto c2 = 0.5411616168555691;
      const auto c3 = 0.2003428171932657;
      const auto c4 = 0.06853891945200943;
      const auto c5 = (137.0/60.0 - clog(-u))/120.0;
      const auto c6 = -1.0/1440.0;

      const double cs[5] = {
         -1.653439153439153e-05, 2.296443268665491e-08,
         -9.941312851365761e-11, 6.691268265342339e-13,
         -5.793305857439255e-15
      };

      return
         cadd(c0,
         cadd(cmul(u, c1),
         cmul(u2, cadd(c2,
         cadd(cmul(u, c3),
         cmul(u2, cadd(c4,
         cadd(cmul(u, c5),
         cmul(u2, cadd(c6,
         cmul(u,  cadd(cs[0],
         cmul(u2, cadd(cs[1],
         cmul(u2, cadd(cs[2],
         cmul(u2, cadd(cs[3],
         cmul(u2, cs[4])))))))))))))))))));
   }

   std::complex<double> u(0.0, 0.0), r(0.0, 0.0);
   double sgn = 1;

   if (az <= 1.0) {
      u = -clog(1.0 - z);
   } else { // az > 1
      const auto arg = pz > 0.0 ? pz - PI : pz + PI;
      const auto lmz = std::complex<double>(lnz, arg); // clog(-z)
      const auto lmz2 = sqr(lmz);
      u = -clog(1.0 - 1.0/z);
      r = -31.0*PI6/15120.0 + lmz2*(-7.0/720.0*PI4 + lmz2*(-1.0/144.0*PI2 - 1.0/720.0*lmz2));
      sgn = -1;
   }

   const auto sum =
      cmul(u, cadd(bf[0],
      cmul(u, cadd(bf[1],
      cmul(u, cadd(bf[2],
      cmul(u, cadd(bf[3],
      cmul(u, cadd(bf[4],
      cmul(u, cadd(bf[5],
      cmul(u, cadd(bf[6],
      cmul(u, cadd(bf[7],
      cmul(u, cadd(bf[8],
      cmul(u, cadd(bf[9],
      cmul(u, cadd(bf[10],
      cmul(u, cadd(bf[11],
      cmul(u, cadd(bf[12],
      cmul(u, cadd(bf[13],
      cmul(u, cadd(bf[14],
      cmul(u, cadd(bf[15],
      cmul(u, cadd(bf[16],
      cmul(u, bf[17])))))))))))))))))))))))))))))))))));

   return sgn*sum + r;
}

/**
 * @brief Complex polylogarithm \f$\mathrm{Li}_6(z)\f$ with long double precision
 * @param z complex argument
 * @return \f$\mathrm{Li}_6(z)\f$
 */
std::complex<long double> Li6(const std::complex<long double>& z)
{
   const long double eps   = 10.0L*std::numeric_limits<long double>::epsilon();
   const long double PI    = 3.14159265358979323846264338327950288L;
   const long double PI2   = PI*PI;
   const long double PI4   = PI2*PI2;
   const long double PI6   = PI2*PI4;
   const long double zeta6 = 1.01734306198444913971451792979092053L;
   const long double bf[] = {
      1.0L,
     -31.0L/64.0L,
      1.52413408779149519890260631001371742e-01L,
     -3.43655558770576131687242798353909465e-02L,
      5.71747972393689986282578875171467764e-03L,
     -6.81804537465706447187928669410150892e-04L,
      4.99603619487344931145170169621327906e-05L,
     -4.91660511960390477202530822164616726e-07L,
     -3.06329751613021637853530304310404212e-07L,
      1.44145992708490953612537421448012923e-08L,
      3.72724382309241065768599245664598825e-09L,
     -3.73008673454876072077977229159261141e-10L,
     -5.12465268160858324340087646115722592e-11L,
      9.05419309566366828868797447620893807e-12L,
      6.73818826155125170776107544938009482e-13L,
     -2.12158311503031353185279689296609530e-13L,
     -6.84088117190116976639204518781909079e-15L,
      4.86911784620055813206387586512917456e-15L,
     -4.84398784998725041650889647473159420e-18L,
     -1.10271048491074909370430302576046237e-16L,
      3.33537969169393816624889411840735964e-18L,
      2.47353074886413529791540987487137046e-18L,
#if LDBL_DIG > 18
     -1.43706164342324920216883687134737009e-19L,
     -5.50471103350981180614826435059791109e-20L,
      4.74677139173272249791309840662617185e-21L,
      1.21583871780681052243739817416294433e-21L,
     -1.41075524035618500414240309078008657e-22L,
     -2.66388312532683465965856437677118103e-23L,
      3.96676574286310079767900226081781935e-24L,
      5.78216973585436153112366193481125963e-25L,
     -1.07877780631642573172998876995922389e-25L,
     -1.24073970867569098990147736977589145e-26L,
      2.87041179178936017042609524684092346e-27L,
      2.62355535630293306165747520383606820e-28L,
     -7.52294854657541272615881214134337672e-29L,
     -5.44017883796246961820291930722306669e-30L,
      1.95025795325101663793862223381656999e-30L,
      1.09784942822051879961178213597012971e-31L,
     -5.01495835741630092074469585415763612e-32L,
     -2.12867375043927610535633806917391780e-33L,
      1.28159440165221259409319852281486752e-33L,
      3.87108447330479441622204568697607460e-35L,
     -3.25941253155837592741689642881678163e-35L,
     -6.25269198847740581093233860701356903e-37L,
      8.25794162051839801918004563317046685e-37L
#endif
   };

   if (is_close(z, 0.0L, eps)) {
      return { 0.0L, 0.0L };
   }
   if (is_close(z, 1.0L, eps)) {
      return { zeta6, 0.0L };
   }
   if (is_close(z, -1.0L, eps)) {
      return { -31.0L*zeta6/32.0L, 0.0L };
   }

   const auto az  = std::abs(z);
   const auto pz  = std::arg(z);
   const auto lnz = std::log(az);

   if (sqr(lnz) + sqr(pz) < 1.0L) { // |log(z)| < 1
      const auto u  = std::complex<long double>(lnz, pz); // clog(z)
      const auto u2 = u*u;
      const auto c0 = zeta6;
      const auto c1 = 1.03692775514336992633136548645703417L; // zeta(5)
      const auto c2 = 0.541161616855569095758001848270583951L;
      const auto c3 = 0.200342817193265714233289693585241665L;
      const auto c4 = 0.0685389194520094348530172986102510496L;
      const auto c5 = (137.0L/60.0L - clog(-u))/120.0L;
      const auto c6 = -1.0L/1440.0L;

      const long double cs[] = {
        -1.65343915343915343915343915343915344e-05L,
         2.29644326866549088771310993533215755e-08L,
        -9.94131285136576141867147158152449158e-11L,
         6.69126826534233941641349048756456164e-13L,
        -5.79330585743925490598570604983944731e-15L,
         5.93014945895224312384148778345583373e-17L,
#if LDBL_DIG > 18
        -6.85052937218694143079665103072690062e-19L,
         8.67589801792722939607545055256974774e-21L,
        -1.18132150428192854833283051865852971e-22L,
         1.70561884258584436997421138705840854e-24L,
        -2.58484268003343989655128609060798194e-26L,
         4.08008103922306292972099603527323696e-28L,
        -6.66771970595289681764999062443512889e-30L,
         1.12276996725126418754155893790528500e-31L,
        -1.94061827643615870372790552830090522e-33L,
         3.43209344566599308274080635472598888e-35L,
        -6.19462586636097236736900573587284992e-37L
#endif
      };

      std::complex<long double> sum(0.0L, 0.0L);

      for (int i = sizeof(cs)/sizeof(cs[0]) - 1; i >= 1; i--) {
         sum = cmul(u2, cadd(cs[i], sum));
      }

      // lowest order terms w/ different powers
      sum = cadd(c0,
         cadd(cmul(u, c1),
         cmul(u2, cadd(c2,
         cadd(cmul(u, c3),
         cmul(u2, cadd(c4,
         cadd(cmul(u, c5),
         cmul(u2, cadd(c6,
         cmul(u,  cadd(cs[0], sum))))))))))));

      return sum;
   }

   std::complex<long double> u(0.0L, 0.0L), r(0.0L, 0.0L);
   long double sgn = 1;

   if (az <= 1.0L) {
      u = -clog(1.0L - z);
   } else { // az > 1
      const auto arg = pz > 0.0 ? pz - PI : pz + PI;
      const auto lmz = std::complex<long double>(lnz, arg); // clog(-z)
      const auto lmz2 = sqr(lmz);
      u = -clog(1.0L - 1.0L/z);
      r = -31.0L*PI6/15120.0L
          + lmz2*(-7.0L/720.0L*PI4 +
                  lmz2*(-1.0L/144.0L*PI2 - 1.0L/720.0L*lmz2));
      sgn = -1;
   }

   std::complex<long double> sum(0.0L, 0.0L);

   for (int i = sizeof(bf)/sizeof(bf[0]) - 1; i >= 0; i--) {
      sum = cmul(u, cadd(bf[i], sum));
   }

   return sgn*sum + r;
}

} // namespace polylogarithm
