// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the MIT License.
// ====================================================================

#include "Li6.hpp"
#include "complex.hpp"
#include <cfloat>
#include <cmath>

namespace polylogarithm {

namespace {

   template <typename T, int N>
   Complex<T> horner(const Complex<T>& z, const T (&coeffs)[N]) noexcept
   {
      static_assert(N >= 2, "more than two coefficients required");

      const T r = z.re + z.re;
      const T s = z.re * z.re + z.im * z.im;
      T a = coeffs[N - 1], b = coeffs[N - 2];

      for (int i = N - 3; i >= 0; --i) {
         const T t = a;
         a = b + r * a;
         b = coeffs[i] - s * t;
      }

      return Complex<T>(z.re*a + b, z.im*a);
   }

} // anonymous namespace

/**
 * @brief Complex polylogarithm \f$\operatorname{Li}_6(z)\f$
 * @param z_ complex argument
 * @return \f$\operatorname{Li}_6(z)\f$
 * @author Alexander Voigt
 */
std::complex<double> Li6(const std::complex<double>& z_) noexcept
{
   const double PI    = 3.1415926535897932;
   const double PI2   = PI*PI;
   const double PI4   = PI2*PI2;
   const double PI6   = PI2*PI4;
   const double zeta6 = 1.0173430619844491;
   const double bf[18] = {
      1.0                   , -31.0/64.0             ,
      1.5241340877914952e-01, -3.4365555877057613e-02,
      5.7174797239368999e-03, -6.8180453746570645e-04,
      4.9960361948734493e-05, -4.9166051196039048e-07,
     -3.0632975161302164e-07,  1.4414599270849095e-08,
      3.7272438230924107e-09, -3.7300867345487607e-10,
     -5.1246526816085832e-11,  9.0541930956636683e-12,
      6.7381882615512517e-13, -2.1215831150303135e-13,
     -6.8408811719011698e-15,  4.8691178462005581e-15
   };

   const Complex<double> z = { std::real(z_), std::imag(z_) };

   if (z.im == 0) {
      if (z.re == 0) {
         return 0.0;
      }
      if (z.re == 1) {
         return zeta6;
      }
      if (z.re == -1) {
         return -31.0*zeta6/32.0;
      }
   }

   const double nz  = norm_sqr(z);
   const double pz  = arg(z);
   const double lnz = 0.5*std::log(nz);

   if (lnz*lnz + pz*pz < 1) { // |log(z)| < 1
      const Complex<double> u(lnz, pz); // log(z)
      const Complex<double> u2 = u*u;
      const double c0 = zeta6;
      const double c1 = 1.0369277551433699; // zeta(5)
      const double c2 = 0.54116161685556910;
      const double c3 = 0.20034281719326571;
      const double c4 = 0.068538919452009435;
      const Complex<double> c5 = (137.0/60.0 - log(-u))/120.0;
      const double c6 = -1.0/1440.0;

      const double cs[5] = {
         -1.6534391534391534e-05, 2.2964432686654909e-08,
         -9.9413128513657614e-11, 6.6912682653423394e-13,
         -5.7933058574392549e-15
      };

      return c0 + u * c1 +
         u2 * (c2 + u * c3 +
         u2 * (c4 + u * c5 +
         u2 * (c6 +
         u * horner(u2, cs))));
   }

   Complex<double> u(0.0, 0.0), rest(0.0, 0.0);
   double sgn = 1;

   if (nz <= 1) {
      u = -log(1.0 - z);
   } else { // nz > 1
      const double arg = pz > 0.0 ? pz - PI : pz + PI;
      const Complex<double> lmz(lnz, arg); // log(-z)
      const Complex<double> lmz2 = lmz*lmz;
      u = -log(1.0 - 1.0/z);
      rest = -31.0*PI6/15120.0 + lmz2*(-7.0/720.0*PI4 + lmz2*(-1.0/144.0*PI2 - 1.0/720.0*lmz2));
      sgn = -1;
   }

   const Complex<double> u2 = u*u;
   const Complex<double> u4 = u2*u2;
   const Complex<double> u8 = u4*u4;

   return
      rest + sgn * (
         u*bf[0] +
         u2*(bf[1] + u*bf[2]) +
         u4*(bf[3] + u*bf[4] + u2*(bf[5] + u*bf[6])) +
         u8*(bf[7] + u*bf[8] + u2*(bf[9] + u*bf[10]) +
             u4*(bf[11] + u*bf[12] + u2*(bf[13] + u*bf[14]))) +
         u8*u8*(bf[15] + u*bf[16] + u2*bf[17])
      );
}

/**
 * @brief Complex polylogarithm \f$\operatorname{Li}_6(z)\f$ with long double precision
 * @param z_ complex argument
 * @return \f$\operatorname{Li}_6(z)\f$
 * @author Alexander Voigt
 */
std::complex<long double> Li6(const std::complex<long double>& z_) noexcept
{
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

   const Complex<long double> z = { std::real(z_), std::imag(z_) };

   if (z.im == 0) {
      if (z.re == 0) {
         return 0.0L;
      }
      if (z.re == 1) {
         return zeta6;
      }
      if (z.re == -1) {
         return -31.0L*zeta6/32.0L;
      }
   }

   const long double nz  = norm_sqr(z);
   const long double pz  = arg(z);
   const long double lnz = 0.5L*std::log(nz);

   if (lnz*lnz + pz*pz < 1) { // |log(z)| < 1
      const Complex<long double> u(lnz, pz); // log(z)
      const Complex<long double> u2 = u*u;
      const long double c0 = zeta6;
      const long double c1 = 1.03692775514336992633136548645703417L; // zeta(5)
      const long double c2 = 0.541161616855569095758001848270583951L;
      const long double c3 = 0.200342817193265714233289693585241665L;
      const long double c4 = 0.0685389194520094348530172986102510496L;
      const Complex<long double> c5 = (137.0L/60.0L - log(-u))/120.0L;
      const long double c6 = -1.0L/1440.0L;

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

      return c0 + u * c1 +
         u2 * (c2 + u * c3 +
         u2 * (c4 + u * c5 +
         u2 * (c6 +
         u * horner(u2, cs))));
   }

   Complex<long double> u(0.0L, 0.0L), rest(0.0L, 0.0L);
   long double sgn = 1;

   if (nz <= 1) {
      u = -log(1.0L - z);
   } else { // nz > 1
      const long double arg = pz > 0.0 ? pz - PI : pz + PI;
      const Complex<long double> lmz(lnz, arg); // log(-z)
      const Complex<long double> lmz2 = lmz*lmz;
      u = -log(1.0L - 1.0L/z);
      rest = -31.0L*PI6/15120.0L
             + lmz2*(-7.0L/720.0L*PI4 +
                     lmz2*(-1.0L/144.0L*PI2 - 1.0L/720.0L*lmz2));
      sgn = -1;
   }

   return rest + sgn*u*horner(u, bf);
}

} // namespace polylogarithm
