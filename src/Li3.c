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


static _Bool is_close(double _Complex a, double b, double eps)
{
   return fabs(creal(a) - b) < eps && fabs(cimag(a)) < eps;
}


static _Bool is_closel(long double _Complex a, long double b, long double eps)
{
   return fabsl(creall(a) - b) < eps && fabsl(cimagl(a)) < eps;
}


/**
 * @brief Complex trilogarithm \f$\mathrm{Li}_3(z)\f$
 * @param z complex argument
 * @return \f$\mathrm{Li}_3(z)\f$
 */
double _Complex cli3(const double _Complex z)
{
   const double eps   = 10.0*DBL_EPSILON;
   const double PI    = 3.141592653589793;
   const double PI2   = PI*PI;
   const double zeta2 = 1.644934066848226;
   const double zeta3 = 1.202056903159594;
   const double bf[18] = {
      1.0                  , -3.0/8.0              ,
      17.0/216.0           , -5.0/576.0            ,
      1.296296296296296e-04,  8.101851851851851e-05,
     -3.419357160853759e-06, -1.328656462585034e-06,
      8.660871756109851e-08,  2.526087595532039e-08,
     -2.144694468364064e-09, -5.140110622012978e-10,
      5.249582114600829e-11,  1.088775440663631e-11,
     -1.277939609449369e-12, -2.369824177308745e-13,
      3.104357887965462e-14,  5.261758629912506e-15
   };

   if (is_close(z, 0.0, eps)) {
      return 0.0;
   }
   if (is_close(z, 1.0, eps)) {
      return zeta3;
   }
   if (is_close(z, -1.0, eps)) {
      return -0.75*zeta3;
   }
   if (is_close(z, 0.5, eps)) {
      const double ln2  = 0.6931471805599453; // ln(2)
      const double ln23 = 0.3330246519889295; // ln(2)^3
      return (-2*PI2*ln2 + 4*ln23 + 21*zeta3)/24.0;
   }

   const double rz  = creal(z);
   const double iz  = cimag(z);
   const double nz  = rz*rz + iz*iz;
   const double pz  = atan2(iz, rz);
   const double lnz = 0.5*log(nz);

   if (lnz*lnz + pz*pz < 1.0) { // |log(z)| < 1
      const double _Complex u  = lnz + pz*I; // clog(z)
      const double _Complex u2 = u*u;
      const double _Complex c0 = zeta3 + u*(zeta2 - u2/12.0);
      const double _Complex c1 = 0.25 * (3.0 - 2.0*fast_clog(-u));

      const double cs[7] = {
         -3.472222222222222e-03, 1.157407407407407e-05,
         -9.841899722852104e-08, 1.148221634332745e-09,
         -1.581572499080917e-11, 2.419500979252515e-13,
         -3.982897776989488e-15
      };

      return c0 +
         u2 * (c1 +
         u2 * (cs[0] +
         u2 * (cs[1] +
         u2 * (cs[2] +
         u2 * (cs[3] +
         u2 * (cs[4] +
         u2 * (cs[5] +
         u2 * (cs[6]))))))));
   }

   double _Complex u = 0.0, rest = 0.0;

   if (nz <= 1.0) {
      u = -fast_clog(1.0 - z);
   } else { // nz > 1
      const double arg = pz > 0.0 ? pz - PI : pz + PI;
      const double _Complex lmz = lnz + arg*I; // clog(-z)
      u = -fast_clog(1.0 - 1.0/z);
      rest = -lmz*(lmz*lmz/6.0 + zeta2);
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
      u * (bf[17]))))))))))))))))));
}


/**
 * @brief Complex trilogarithm \f$\mathrm{Li}_3(z)\f$
 * @param z complex argument
 * @return \f$\mathrm{Li}_3(z)\f$
 */
long double _Complex cli3l(const long double _Complex z)
{
   const long double eps   = 10.0L*LDBL_EPSILON;
   const long double PI    = 3.14159265358979323846264338327950288L;
   const long double PI2   = PI*PI;
   const long double zeta2 = 1.64493406684822643647241516664602519L;
   const long double zeta3 = 1.20205690315959428539973816151144999L;
   const long double bf[] = {
      1.0L,
     -3.0L/8.0L,
      17.0L/216.0L,
     -5.0L/576.0L,
      7.0L/54000.0L,
      7.0L/86400.0L,
     -3.41935716085375949321527552820069827e-06L,
     -1.32865646258503401360544217687074830e-06L,
      8.66087175610985134794658604182413706e-08L,
      2.52608759553203997648442092886537331e-08L,
     -2.14469446836406476093388507573649032e-09L,
     -5.14011062201297891533581769272004962e-10L,
      5.24958211460082943639408880855807284e-11L,
      1.08877544066363183753729715704249107e-11L,
     -1.27793960944936953055818317540722120e-12L,
     -2.36982417730874520997977788101244891e-13L,
      3.10435788796546229428475327046556211e-14L,
      5.26175862991250608413183925112250061e-15L,
     -7.53847954994926536599250143226771028e-16L,
     -1.18623225777522852530825009512459322e-16L,
      1.83169799654913833820892731212815349e-17L,
      2.70681710318373501514907347126169436e-18L,
#if LDBL_DIG > 18
     -4.45543389782963882643263099217632212e-19L,
     -6.23754849225569465036532224739838641e-20L,
      1.08515215348745349131365609968642833e-20L,
      1.44911748660360819307349049665275324e-21L,
     -2.64663397544589903347408911861443741e-22L,
     -3.38976534885101047219258165860814078e-23L,
      6.46404773360331088903253098219534234e-24L,
      7.97583448960241242420922272590502795e-25L,
     -1.58091787902874833559211176293826770e-25L,
     -1.88614997296228681931102253988531956e-26L,
      3.87155366384184733039971271888313319e-27L,
      4.48011750023456073048653898320511684e-28L,
     -9.49303387191183612641753676027699150e-29L,
     -1.06828138090773812240182143033807908e-29L,
      2.33044789361030518600785199019281371e-30L,
      2.55607757265197540805635698286695865e-31L,
     -5.72742160613725968447274458033057100e-32L,
     -6.13471321379642358258549296897773326e-33L,
      1.40908086040689448401268688489421700e-33L,
      1.47642223976665341443182801167106626e-34L,
     -3.47010516489959160555004020312910903e-35L,
     -3.56210662409746357967357370318293608e-36L,
      8.55369656823692105754731289124468101e-37L
#endif
   };

   if (is_closel(z, 0.0L, eps)) {
      return 0.0L;
   }
   if (is_closel(z, 1.0L, eps)) {
      return zeta3;
   }
   if (is_closel(z, -1.0L, eps)) {
      return -0.75L*zeta3;
   }
   if (is_closel(z, 0.5L, eps)) {
      const long double ln2  = 0.693147180559945309417232121458176568L; // ln(2)
      const long double ln23 = 0.333024651988929479718853582611730544L; // ln(2)^3
      return (-2*PI2*ln2 + 4*ln23 + 21*zeta3)/24.0L;
   }

   const long double rz  = creall(z);
   const long double iz  = cimagl(z);
   const long double nz  = rz*rz + iz*iz;
   const long double pz  = atan2l(iz, rz);
   const long double lnz = 0.5L*logl(nz);

   if (lnz*lnz + pz*pz < 1.0L) { // |log(z)| < 1
      const long double _Complex u  = lnz + pz*I; // clog(z)
      const long double _Complex u2 = u*u;
      const long double _Complex c0 = zeta3 + u*(zeta2 - u2/12.0L);
      const long double _Complex c1 = 0.25L * (3.0L - 2.0L*fast_clogl(-u));

      const long double cs[] = {
        -3.47222222222222222222222222222222222e-03L,
         1.15740740740740740740740740740740741e-05L,
        -9.84189972285210380448475686570924666e-08L,
         1.14822163433274544385655496766607878e-09L,
        -1.58157249908091658933409775160616911e-11L,
         2.41950097925251519452732701564998016e-13L,
        -3.98289777698948774786517290926462002e-15L,
         6.92336661830592905806820954095065870e-17L,
        -1.25527223044997727545846570912655367e-18L,
#if LDBL_DIG > 18
         2.35375400276846523056441171414060379e-20L,
        -4.53639890345868701844750708901700830e-22L,
         8.94516967039264316712031170773304472e-24L,
        -1.79828400469549627172020247141015426e-25L,
         3.67549976479373844433604733912674099e-27L,
        -7.62080797156479522953948500963765478e-29L,
         1.60004196436948597517376392257325602e-30L,
        -3.39676114756037558792312060520851852e-32L,
         7.28227228675776469531725636144432664e-34L,
        -1.57502264795800348718497893940378261e-35L,
         3.43354009248058933588797212016527038e-37L
#endif
      };

      long double _Complex sum = 0.0L;

      for (int i = sizeof(cs)/sizeof(cs[0]) - 1; i >= 0; i--) {
         sum = u2 * (cs[i] + sum);
      }

      // lowest order terms w/ different powers
      sum = c0 + u2 * (c1 + sum);

      return sum;
   }

   long double _Complex u = 0.0L, rest = 0.0L;

   if (nz <= 1.0L) {
      u = -fast_clogl(1.0L - z);
   } else { // nz > 1
      const long double arg = pz > 0.0L ? pz - PI : pz + PI;
      const long double _Complex lmz = lnz + arg*I; // clog(-z)
      u = -fast_clogl(1.0L - 1.0L/z);
      rest = -lmz*(lmz*lmz/6.0L + zeta2);
   }

   long double _Complex sum = 0.0L;

   for (int i = sizeof(bf)/sizeof(bf[0]) - 1; i >= 0; i--) {
      sum = u * (bf[i] + sum);
   }

   return rest + sum;
}


/** C++ wrapper for cli2 */
void cli3_(double re, double im, double* res_re, double* res_im)
{
   const double _Complex result = cli3(re + I*im);
   *res_re = creal(result);
   *res_im = cimag(result);
}


/** C++ wrapper for cli2 */
void cli3l_(long double re, long double im, long double* res_re, long double* res_im)
{
   const long double _Complex result = cli3l(re + I*im);
   *res_re = creall(result);
   *res_im = cimagl(result);
}
