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


static _Bool is_close(double _Complex a, double b, double eps)
{
   return fabs(creal(a) - b) < eps && fabs(cimag(a)) < eps;
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
      return 0.0 + 0.0*I;
   }
   if (is_close(z, 1.0, eps)) {
      return zeta3 + 0.0*I;
   }
   if (is_close(z, -1.0, eps)) {
      return -0.75*zeta3 + 0.0*I;
   }
   if (is_close(z, 0.5, eps)) {
      const double ln2  = 0.6931471805599453; // ln(2)
      const double ln23 = 0.3330246519889295; // ln(2)^3
      return (-2*PI2*ln2 + 4*ln23 + 21*zeta3)/24.0 + 0.0*I;
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

   double _Complex u = 0.0 + 0.0*I, rest = 0.0 + 0.0*I;

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


/** C++ wrapper for cli2 */
void cli3_(double re, double im, double* res_re, double* res_im)
{
   const double _Complex result = cli3(re + I*im);
   *res_re = creal(result);
   *res_im = cimag(result);
}
