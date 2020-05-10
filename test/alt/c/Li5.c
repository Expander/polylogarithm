/* ====================================================================
 * This file is part of Polylogarithm.
 *
 * Polylogarithm is licenced under the GNU Lesser General Public
 * License (GNU LGPL) version 3.
 * ==================================================================== */

#include <complex.h>
#include <float.h>
#include <math.h>


static double _Complex fast_clog(double _Complex z)
{
   const double rz = creal(z) == 0.0 ? fabs(creal(z)) : creal(z);
   const double iz = cimag(z) == 0.0 ? fabs(cimag(z)) : cimag(z);

   return 0.5*log(rz*rz + iz*iz) + I*atan2(iz, rz);
}


static long double _Complex fast_clogl(long double _Complex z)
{
   const long double rz = creall(z) == 0.0L ? fabsl(creall(z)) : creall(z);
   const long double iz = cimagl(z) == 0.0L ? fabsl(cimagl(z)) : cimagl(z);

   return 0.5L*logl(rz*rz + iz*iz) + I*atan2l(iz, rz);
}


/**
 * @brief Complex polylogarithm \f$\mathrm{Li}_5(z)\f$
 * @param z complex argument
 * @return \f$\mathrm{Li}_5(z)\f$
 */
double _Complex cli5(double _Complex z)
{
   const double PI    = 3.141592653589793;
   const double PI2   = PI*PI;
   const double PI4   = PI2*PI2;
   const double zeta5 = 1.036927755143369;
   const double bf[19] = {
      1.0, -15.0/32.0,
      1.395318930041152e-01, -2.863377700617283e-02,
      4.031741255144032e-03, -3.398501800411522e-04,
      4.544518462161766e-06,  2.391680804856901e-06,
     -1.276269260012274e-07, -3.162898430650593e-08,
      3.284811844533519e-09,  4.761371399566057e-10,
     -8.084689817190983e-11, -7.238764858773720e-12,
      1.943976011517396e-12,  1.025697840597723e-13,
     -4.618055100988483e-14, -1.153585719647058e-15,
      1.090354540133339e-15
   };

   const double rz  = creal(z);
   const double iz  = cimag(z);

   if (iz == 0) {
      if (rz == 0) {
         return 0.0;
      }
      if (rz == 1) {
         return zeta5;
      }
      if (rz == -1) {
         return -15.0*zeta5/16.0;
      }
   }

   const double nz  = rz*rz + iz*iz;
   const double pz  = atan2(iz, rz);
   const double lnz = 0.5*log(nz);

   if (lnz*lnz + pz*pz < 1.0) { // |log(z)| < 1
      const double _Complex u  = lnz + pz*I; // clog(z)
      const double _Complex u2 = u*u;
      const double c0 = zeta5;
      const double c1 = 1.082323233711138; // zeta(4)
      const double c2 = 0.6010284515797971; // zeta(3)/2
      const double c3 = 0.2741556778080377;
      const double _Complex c4 = (25.0/12.0 - fast_clog(-u))/24.0;
      const double c5 = -1.0/240.0;

      const double cs[6] = {
         -1.157407407407407e-04, 2.066798941798942e-07,
         -1.093544413650234e-09, 8.698648744945041e-12,
         -8.689958786158882e-14, 1.008125408021881e-15
      };

      return c0 + u * c1 +
         u2 * (c2 + u * c3 +
         u2 * (c4 + u * c5 +
         u2 * (cs[0] +
         u2 * (cs[1] +
         u2 * (cs[2] +
         u2 * (cs[3] +
         u2 * (cs[4] +
         u2 * (cs[5]))))))));
   }

   double _Complex u = 0.0, rest = 0.0;

   if (nz <= 1.0) {
      u = -fast_clog(1.0 - z);
   } else { // nz > 1
      const double arg = pz > 0.0 ? pz - PI : pz + PI;
      const double _Complex lmz = lnz + arg*I; // clog(-z)
      const double _Complex lmz2 = lmz*lmz;
      u = -fast_clog(1.0 - 1.0/z);
      rest = -1.0/360.0*lmz*(7*PI4 + lmz2*(10.0*PI2 + 3.0*lmz2));
   }

   return rest +
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
      u * (bf[18])))))))))))))))))));
}


/**
 * @brief Complex polylogarithm \f$\mathrm{Li}_5(z)\f$ with long double precision
 * @param z complex argument
 * @return \f$\mathrm{Li}_5(z)\f$
 */
long double _Complex cli5l(long double _Complex z)
{
   const long double PI    = 3.14159265358979323846264338327950288L;
   const long double PI2   = PI*PI;
   const long double PI4   = PI2*PI2;
   const long double zeta5 = 1.03692775514336992633136548645703417L;
   const long double bf[] = {
      1.0L,
     -15.0L/32.0L,
      1.39531893004115226337448559670781893e-01L,
     -2.86337770061728395061728395061728395e-02L,
      4.03174125514403292181069958847736626e-03L,
     -3.39850180041152263374485596707818930e-04L,
      4.54451846216176664909446003743133026e-06L,
      2.39168080485690118829088702752453967e-06L,
     -1.27626926001227465885443518981747128e-07L,
     -3.16289843065059324402567872795470007e-08L,
      3.28481184453351916215185742384719818e-09L,
      4.76137139956605790483191328977324967e-10L,
     -8.08468981719098302564602603623317490e-11L,
     -7.23876485877372069468292158375897580e-12L,
      1.94397601151739684930556492093599766e-12L,
      1.02569784059772359718813559433198981e-13L,
     -4.61805510098848301805820862410656158e-14L,
     -1.15358571964705800368425114834190972e-15L,
      1.09035454013333939879770883809662643e-15L,
      2.31481363172925263940797103190091493e-18L,
     -2.56699170432652921943348919933966693e-17L,
#if LDBL_DIG > 18
      4.57086206073149690144959626860139115e-19L,
      6.03667796132057058823561033114107090e-19L,
     -2.16776249440624129587941717218396578e-20L,
     -1.41940966156001652983322668820112130e-20L,
      7.50200095064138625532377521619527234e-22L,
      3.33870453950783971643715159254469304e-22L,
     -2.30600404426203476825215151352586388e-23L,
     -7.85817324568948189044990646315027350e-24L,
      6.66834530437388085486513704613056895e-25L,
      1.85091565409252971894649796883651603e-25L,
     -1.85915294451740855841031840576364891e-26L,
     -4.36297464803458904472660817437095794e-27L,
      5.06110760995292844822634895878109349e-28L,
      1.02919182497568782037888979300008731e-28L,
     -1.35513912210183166156877853283041765e-29L,
     -2.42940596129573826559241540956570701e-30L,
      3.58519739665037052115164738053066333e-31L,
      5.73796581610397206400846538673280837e-32L,
     -9.40035936245687345352774520389480989e-33L,
     -1.35590280493486311500090284171223034e-33L,
      2.44784384191528918377859141748058708e-34L,
      3.20528849130720958124170712217928968e-35L,
     -6.33983878185254827964152375942174520e-36L,
     -7.57925545801218291534870941639851689e-37L
#endif
   };

   const long double rz  = creall(z);
   const long double iz  = cimagl(z);

   if (iz == 0) {
      if (rz == 0) {
         return 0.0L;
      }
      if (rz == 1) {
         return zeta5;
      }
      if (rz == -1) {
         return -15.0L*zeta5/16.0L;
      }
   }

   const long double nz  = rz*rz + iz*iz;
   const long double pz  = atan2l(iz, rz);
   const long double lnz = 0.5L*logl(nz);

   if (lnz*lnz + pz*pz < 1.0L) { // |log(z)| < 1
      const long double _Complex u  = lnz + pz*I; // clog(z)
      const long double _Complex u2 = u*u;
      const long double c0 = zeta5;
      const long double c1 = 1.08232323371113819151600369654116790L; // zeta(4)
      const long double c2 = 0.601028451579797142699869080755724995L; // zeta(3)/2
      const long double c3 = 0.274155677808037739412069194441004198L;
      const long double _Complex c4 = (25.0L/12.0L - fast_clogl(-u))/24.0L;
      const long double c5 = -1.0L/240.0L;

      const long double cs[] = {
        -1.15740740740740740740740740740740741e-04L,
         2.06679894179894179894179894179894180e-07L,
        -1.09354441365023375605386187396769407e-09L,
         8.69864874494504124133753763383393013e-12L,
        -8.68995878615888235897855907475917096e-14L,
         1.00812540802188133105305292318749173e-15L,
        -1.30160058071551887185136369583811112e-17L,
#if LDBL_DIG > 18
         1.82193858376471817317584461603964703e-19L,
        -2.71703945984843566116551019291461834e-21L,
         4.26404710646461092493552846764602135e-23L,
        -6.97907523609028772068847244464155123e-25L,
         1.18322350137468824961908885022923872e-26L,
        -2.06699310884539801347149709357488995e-28L,
         3.70514089192917181888714449508744051e-30L,
        -6.79216396752655546304766934905316826e-32L,
         1.26987457489641744061409835124861589e-33L,
        -2.41590408788077922327391223699041147e-35L,
         4.66812326074215685597260023169508118e-37L
#endif
      };

      long double _Complex sum = 0.0L;

      for (int i = sizeof(cs)/sizeof(cs[0]) - 1; i >= 0; i--) {
         sum = u2 *(cs[i] + sum);
      }

      return c0 + u * c1 +
         u2 * (c2 + u * c3 +
         u2 * (c4 + u * c5 + sum));
   }

   long double _Complex u = 0.0L, rest = 0.0L;

   if (nz <= 1.0L) {
      u = -fast_clogl(1.0L - z);
   } else { // nz > 1
      const long double arg = pz > 0.0 ? pz - PI : pz + PI;
      const long double _Complex lmz = lnz + arg*I; // clog(-z)
      const long double _Complex lmz2 = lmz*lmz;
      u = -fast_clogl(1.0L - 1.0L/z);
      rest = -1.0L/360.0L*lmz*(7*PI4 + lmz2*(10.0L*PI2 + 3.0L*lmz2));
   }

   long double _Complex sum = 0.0L;

   for (int i = sizeof(bf)/sizeof(bf[0]) - 1; i >= 0; i--) {
      sum = u * (bf[i] + sum);
   }

   return rest + sum;
}


/** C++ wrapper for cli5 */
void cli5_(double re, double im, double* res_re, double* res_im)
{
   const double _Complex result = cli5(re + I*im);
   *res_re = creal(result);
   *res_im = cimag(result);
}


/** C++ wrapper for cli5l */
void cli5l_(long double re, long double im, long double* res_re, long double* res_im)
{
   const long double _Complex result = cli5l(re + I*im);
   *res_re = creall(result);
   *res_im = cimagl(result);
}
