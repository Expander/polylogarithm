// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the GNU Lesser General Public
// License (GNU LGPL) version 3.
// ====================================================================

#include "Li2.hpp"
#include <cfloat>
#include <limits>

namespace polylogarithm {

namespace {
   template <typename T>
   std::complex<T> clog(const std::complex<T>& z) noexcept
   {
      const T rz = std::real(z);
      const T iz = std::imag(z);

      return std::complex<T>(0.5*std::log(rz*rz + iz*iz), std::atan2(iz, rz));
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

   template <int Nstart, int Nend, typename T, int N>
   std::complex<T> horner(const std::complex<T>& z, const T (&coeffs)[N]) noexcept
   {
      static_assert(Nstart <= Nend && Nend < N && Nend >= 1, "invalid array bounds");

      const T x = std::real(z);
      const T y = std::imag(z);
      const T r = x + x;
      const T s = x * x + y * y;
      T a = coeffs[Nend], b = coeffs[Nend - 1];

      for (int i = Nend - 2; i >= Nstart; --i) {
         const T t = a;
         a = b + r * a;
         b = coeffs[i] - s * t;
      }

      return std::complex<T>(x*a + b, y*a);
   }

} // anonymous namespace

/**
 * @brief Real dilogarithm \f$\mathrm{Li}_2(x)\f$
 * @param x real argument
 * @return \f$\mathrm{Li}_2(x)\f$
 * @note Implementation translated by R.Brun from CERNLIB DILOG function C332
 * @author K.S. Kölbig
 *
 * Implemented as a truncated series expansion in terms of Chebyshev
 * polynomials, see [Yudell L. Luke: Mathematical functions and their
 * approximations, Academic Press Inc., New York 1975, p.67].
 */
double Li2(double x) noexcept
{
   const double PI  = 3.141592653589793;
   const double HF  = 0.5;
   const double PI2 = PI*PI;
   const double PI3 = PI2/3;
   const double PI6 = PI2/6;
   const double PI12 = PI2/12;
   const double C[20] = {  0.42996693560813697, 0.40975987533077106,
     -0.01858843665014592, 0.00145751084062268,-0.00014304184442340,
      0.00001588415541880,-0.00000190784959387, 0.00000024195180854,
     -0.00000003193341274, 0.00000000434545063,-0.00000000060578480,
      0.00000000008612098,-0.00000000001244332, 0.00000000000182256,
     -0.00000000000027007, 0.00000000000004042,-0.00000000000000610,
      0.00000000000000093,-0.00000000000000014, 0.00000000000000002};

   double T = 0, H = 0, Y = 0, S = 0, A = 0, ALFA = 0, B1 = 0, B2 = 0, B0 = 0;

   if (x == 1) {
       H = PI6;
   } else if (x == -1) {
       H = -PI12;
   } else {
       T = -x;
       if (T <= -2) {
           Y = -1/(1+T);
           S = 1;
           B1= std::log(-T);
           B2= std::log(1+1/T);
           A = -PI3+HF*(B1*B1-B2*B2);
       } else if (T < -1) {
           Y = -1-T;
           S = -1;
           A = std::log(-T);
           A = -PI6+A*(A+std::log(1+1/T));
       } else if (T <= -0.5) {
           Y = -(1+T)/T;
           S = 1;
           A = std::log(-T);
           A = -PI6+A*(-HF*A+std::log(1+T));
       } else if (T < 0) {
           Y = -T/(1+T);
           S = -1;
           B1= std::log(1+T);
           A = HF*B1*B1;
       } else if (T <= 1) {
           Y = T;
           S = 1;
           A = 0;
       } else {
           Y = 1/T;
           S = -1;
           B1= std::log(T);
           A = PI6+HF*B1*B1;
       }
       H    = Y+Y-1;
       ALFA = H+H;
       B1   = 0;
       B2   = 0;
       for (int i = 19; i >= 0; i--) {
          B0 = C[i] + ALFA*B1-B2;
          B2 = B1;
          B1 = B0;
       }
       H = -(S*(B0-H*B2)+A);
    }
    return H;
}

/**
 * @brief Real dilogarithm \f$\mathrm{Li}_2(z)\f$ with long double precision
 * @param x real argument
 * @return \f$\mathrm{Li}_2(z)\f$
 * @author K.S. Kölbig
 * @note Implementation based on translation by R.Brun from CERNLIB
 *    DILOG function C332, extended by Alexander Voigt to long double
 *    precision
 *
 * Implemented as a truncated series expansion in terms of Chebyshev
 * polynomials, see [Yudell L. Luke: Mathematical functions and their
 * approximations, Academic Press Inc., New York 1975, p.67].
 */
long double Li2(long double x) noexcept
{
   const long double PI  = 3.14159265358979323846264338327950288L;
   const long double HF  = 0.5L;
   const long double PI2 = PI*PI;
   const long double PI3 = PI2/3;
   const long double PI6 = PI2/6;
   const long double PI12 = PI2/12;
   const long double C[] = {
      0.4299669356081369720370336786993879912L,
      0.4097598753307710584682637109252552781L,
     -0.0185884366501459196476416491402122676L,
      0.0014575108406226785536739284164594927L,
     -0.0001430418444234004877488301200908765L,
      0.0000158841554187955323619055047167740L,
     -0.0000019078495938658272271963211420884L,
      0.0000002419518085416474949946146434290L,
     -0.0000000319334127425178346049601414286L,
      0.0000000043454506267691229879571784782L,
     -0.0000000006057848011840744442970533091L,
      0.0000000000861209779935949824428368452L,
     -0.0000000000124433165993886798964242137L,
      0.0000000000018225569623573633006554774L,
     -0.0000000000002700676604911465180857223L,
      0.0000000000000404220926315266464883286L,
     -0.0000000000000061032514526918795037783L,
      0.0000000000000009286297533019575861303L,
     -0.0000000000000001422602085511244683975L,
      0.0000000000000000219263171815395735398L,
     -0.0000000000000000033979732421589786340L,
#if LDBL_DIG > 18
      0.0000000000000000005291954244833147146L,
     -0.0000000000000000000827858081427899765L,
      0.0000000000000000000130037173454556037L,
     -0.0000000000000000000020502222425528249L,
      0.0000000000000000000003243578549148930L,
     -0.0000000000000000000000514779990334321L,
      0.0000000000000000000000081938774771716L,
     -0.0000000000000000000000013077835405713L,
      0.0000000000000000000000002092562930580L,
     -0.0000000000000000000000000335616615054L,
      0.0000000000000000000000000053946577714L,
     -0.0000000000000000000000000008689193209L,
      0.0000000000000000000000000001402281687L,
     -0.0000000000000000000000000000226715578L,
      0.0000000000000000000000000000036717417L,
     -0.0000000000000000000000000000005956152L,
      0.0000000000000000000000000000000967662L,
     -0.0000000000000000000000000000000157439L,
      0.0000000000000000000000000000000025650L,
     -0.0000000000000000000000000000000004185L,
      0.0000000000000000000000000000000000683L,
     -0.0000000000000000000000000000000000112L,
      0.0000000000000000000000000000000000018L,
     -0.0000000000000000000000000000000000003L
#endif
   };

   long double T = 0, H = 0, Y = 0, S = 0, A = 0, ALFA = 0, B1 = 0, B2 = 0, B0 = 0;

   if (x == 1) {
       H = PI6;
   } else if (x == -1) {
       H = -PI12;
   } else {
       T = -x;
       if (T <= -2) {
           Y = -1/(1+T);
           S = 1;
           B1= std::log(-T);
           B2= std::log(1+1/T);
           A = -PI3+HF*(B1*B1-B2*B2);
       } else if (T < -1) {
           Y = -1-T;
           S = -1;
           A = std::log(-T);
           A = -PI6+A*(A+std::log(1+1/T));
       } else if (T <= -0.5L) {
           Y = -(1+T)/T;
           S = 1;
           A = std::log(-T);
           A = -PI6+A*(-HF*A+std::log(1+T));
       } else if (T < 0) {
           Y = -T/(1+T);
           S = -1;
           B1= std::log(1+T);
           A = HF*B1*B1;
       } else if (T <= 1) {
           Y = T;
           S = 1;
           A = 0;
       } else {
           Y = 1/T;
           S = -1;
           B1= std::log(T);
           A = PI6+HF*B1*B1;
       }
       H    = Y+Y-1;
       ALFA = H+H;
       B1   = 0;
       B2   = 0;
       for (int i = sizeof(C)/sizeof(C[0]) - 1; i >= 0; i--) {
          B0 = C[i] + ALFA*B1-B2;
          B2 = B1;
          B1 = B0;
       }
       H = -(S*(B0-H*B2)+A);
    }
    return H;
}

/**
 * @brief Complex dilogarithm \f$\mathrm{Li}_2(z)\f$
 * @param z complex argument
 * @return \f$\mathrm{Li}_2(z)\f$
 * @note Implementation translated from SPheno to C++
 * @author Werner Porod
 * @note translated to C++ by Alexander Voigt
 */
std::complex<double> Li2(const std::complex<double>& z) noexcept
{
   const double PI = 3.141592653589793;

   // bf[1..N-1] are the even Bernoulli numbers / (2 n + 1)!
   // generated by: Table[BernoulliB[2 n]/(2 n + 1)!, {n, 1, 9}]
   const double bf[10] = {
      - 1.0/4.0,
      + 1.0/36.0,
      - 1.0/3600.0,
      + 1.0/211680.0,
      - 1.0/10886400.0,
      + 1.0/526901760.0,
      - 4.064761645144226e-11,
      + 8.921691020456453e-13,
      - 1.993929586072108e-14,
      + 4.518980029619918e-16
   };

   const double rz = std::real(z);
   const double iz = std::imag(z);
   const double nz = rz*rz + iz*iz;

   // special cases
   if (iz == 0.0) {
      if (rz <= 1.0) {
         return Li2(rz);
      }
      if (rz > 1.0) {
         return std::complex<double>(Li2(rz), -PI*std::log(rz));
      }
   } else if (nz < std::numeric_limits<double>::epsilon()) {
      return z;
   }

   std::complex<double> cy(0.0, 0.0), cz(0.0, 0.0);
   int jsgn = 0, ipi12 = 0;

   // transformation to |z|<1, Re(z)<=0.5
   if (rz <= 0.5) {
      if (nz > 1.0) {
         const std::complex<double> lz = clog(-z);
         cy = -0.5 * lz*lz;
         cz = -clog(1.0 - 1.0 / z);
         jsgn = -1;
         ipi12 = -2;
      } else { // nz <= 1
         cy = 0;
         cz = -clog(1.0 - z);
         jsgn = 1;
         ipi12 = 0;
      }
   } else { // rz > 0.5
      if (nz <= 2*rz) {
         cz = -clog(z);
         cy = cz * clog(1.0 - z);
         jsgn = -1;
         ipi12 = 2;
      } else { // nz > 2*rz
         const std::complex<double> lz = clog(-z);
         cy = -0.5 * lz*lz;
         cz = -clog(1.0 - 1.0 / z);
         jsgn = -1;
         ipi12 = -2;
      }
   }

   // the dilogarithm
   const std::complex<double> cz2(cz*cz);

   std::complex<double> sum = horner<1, 9>(cz2, bf);

   // lowest order terms w/ different powers
   sum = cadd(cz, cmul(cz2, cadd(bf[0], cmul(cz, sum))));

   return double(jsgn) * sum + cy + ipi12 * PI * PI / 12.0;
}

/**
 * @brief Complex dilogarithm \f$\mathrm{Li}_2(z)\f$ with long double precision
 * @param z complex argument
 * @return \f$\mathrm{Li}_2(z)\f$
 * @note Implementation translated from SPheno to C++
 * @author Werner Porod
 * @note translated to C++ and extended to long double precision by Alexander Voigt
 */
std::complex<long double> Li2(const std::complex<long double>& z) noexcept
{
   const long double PI = 3.14159265358979323846264338327950288L;

   // bf[1..N-1] are the even Bernoulli numbers / (2 n + 1)!
   // generated by: Table[BernoulliB[2 n]/(2 n + 1)!, {n, 1, 22}]
   const long double bf[] = {
      -1.0L/4.0L                                 ,
       1.0L/36.0L                                ,
      -1.0L/3600.0L                              ,
       1.0L/211680.0L                            ,
      -1.0L/10886400.0L                          ,
       1.0L/526901760.0L                         ,
      -4.06476164514422552680590938629196667e-11L,
       8.92169102045645255521798731675274885e-13L,
      -1.99392958607210756872364434779378971e-14L,
       4.51898002961991819165047655285559323e-16L,
      -1.03565176121812470144834115422186567e-17L,
#if LDBL_DIG > 18
       2.39521862102618674574028374300098038e-19L,
      -5.58178587432500933628307450562541991e-21L,
       1.30915075541832128581230739918659230e-22L,
      -3.08741980242674029324227976486646243e-24L,
       7.31597565270220342035790560925214859e-26L,
      -1.74084565723400074098905514775970255e-27L,
       4.15763564461389971961789962077522667e-29L,
      -9.96214848828462210319400670245583885e-31L,
       2.39403442489616530052116798789374956e-32L,
      -5.76834735536739008429179316187765424e-34L,
       1.39317947964700797782788660391154833e-35L,
      -3.37212196548508947046847363525493096e-37L
#endif
   };

   constexpr int N = sizeof(bf)/sizeof(bf[0]);

   const long double rz = std::real(z);
   const long double iz = std::imag(z);
   const long double nz = rz*rz + iz*iz;

   // special cases
   if (iz == 0.0L) {
      if (rz <= 1.0L) {
         return Li2(rz);
      }
      if (rz > 1.0L) {
         return std::complex<long double>(Li2(rz), -PI*std::log(rz));
      }
   } else if (nz < std::numeric_limits<long double>::epsilon()) {
      return z;
   }

   std::complex<long double> cy(0.0L, 0.0L), cz(0.0L, 0.0L);
   int jsgn = 0, ipi12 = 0;

   // transformation to |z|<1, Re(z)<=0.5
   if (rz <= 0.5L) {
      if (nz > 1.0L) {
         const std::complex<long double> lz = clog(-z);
         cy = -0.5L * lz*lz;
         cz = -clog(1.0L - 1.0L/z);
         jsgn = -1;
         ipi12 = -2;
      } else { // nz <= 1.0L
         cy = 0;
         cz = -clog(1.0L - z);
         jsgn = 1;
         ipi12 = 0;
      }
   } else { // rz > 0.5L
      if (nz <= 2*rz) {
         cz = -clog(z);
         cy = cz * clog(1.0L - z);
         jsgn = -1;
         ipi12 = 2;
      } else { // nz > 2*rz
         const std::complex<long double> lz = clog(-z);
         cy = -0.5L * lz*lz;
         cz = -clog(1.0L - 1.0L/z);
         jsgn = -1;
         ipi12 = -2;
      }
   }

   const std::complex<long double> cz2(cz*cz);

   std::complex<long double> sum = horner<1, N-1>(cz2, bf);

   // lowest order terms w/ different powers
   sum = cadd(cz, cmul(cz2, cadd(bf[0], cmul(cz, sum))));

   return static_cast<long double>(jsgn)*sum + cy + ipi12*PI*PI/12.0L;
}

} // namespace polylogarithm
