#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <float.h>
#include "tsil.h"

TSIL_REAL TSIL_dilog_real           (TSIL_REAL);
TSIL_REAL TSIL_dilog_series_real    (TSIL_REAL);
TSIL_REAL TSIL_dilog_largeabs_real  (TSIL_REAL);
TSIL_REAL TSIL_dilog_neg_real       (TSIL_REAL);
TSIL_REAL TSIL_dilog_CLZseries_real (TSIL_REAL);

int TSIL_dilog_series_complex    (TSIL_REAL, TSIL_REAL, TSIL_REAL *, TSIL_REAL *);
int TSIL_dilog_largeabs_complex  (TSIL_REAL, TSIL_REAL, TSIL_REAL *, TSIL_REAL *);
int TSIL_dilog_neg_complex       (TSIL_REAL, TSIL_REAL, TSIL_REAL *, TSIL_REAL *);
int TSIL_dilog_CLZseries_complex (TSIL_REAL, TSIL_REAL, TSIL_REAL *, TSIL_REAL *);

#define zeta2L    1.64493406684822643647241516665L
#define M_PI_LONG 3.14159265358979323846264338328L

#define dilog_crossover_high 160.L
#define dilog_crossover_low  0.04L

/* ************************************************************** */

TSIL_REAL TSIL_dilog_real (TSIL_REAL x)
{
  TSIL_REAL xsquared = x * x;

  if (x == 1.0L)
    return zeta2L;
  else if (x == 0.0L)
    return 0.0L;
  else if (xsquared < dilog_crossover_low)
    return TSIL_dilog_series_real (x);
  else if (xsquared > dilog_crossover_high)
    return TSIL_dilog_largeabs_real (x);
  else if (x < 0.0L) 
    return TSIL_dilog_neg_real (x);
  else
    return TSIL_dilog_CLZseries_real (x);
}

/* ************************************************************** */
/* 
   Assumes |x| < 1, and is only called here for 
   |x|^2 < smallcrossover.
*/

TSIL_REAL TSIL_dilog_series_real (TSIL_REAL x)
{
  TSIL_REAL xsquared = x * x;
  TSIL_REAL xcubed, xtothek;
  TSIL_REAL termk;
  TSIL_REAL ksquared;
  TSIL_REAL result = 0.0L;
  int kmax, k;

  xcubed = xsquared * x;
  xtothek = xcubed;
  kmax = ceil( (float) (TSIL_LOG(TSIL_TOL)/TSIL_LOG(TSIL_FABS(x))) ) + 1;
  if (kmax < 5)
    kmax = 5;

  /* Add up the smaller terms first, to help avoid roundoff errors. */
  for (k = 4; k < kmax; k++)
    {
      xtothek *= x;
      ksquared = (TSIL_REAL) (k * k);
      termk = xtothek / ksquared;
      result += termk;
    }

  result += xcubed / 9.0L;
  result += 0.25L * xsquared;
  result += x;

  return result;
}

/* ************************************************************** */
/* 
   The following deals with |x| > 1 by implementing the identity: Li_2
   (x) = -Li_2 (1/x) -zeta(2) - (1/2) (log(-x))^2 It should only be
   called here for |x| > largecrossover.
*/

TSIL_REAL TSIL_dilog_largeabs_real (TSIL_REAL x)
{
  TSIL_REAL logminusx_re;
  TSIL_REAL t1, t2;
  TSIL_REAL xsquared = x * x;
  TSIL_REAL xinv = 1.0L / x;
  TSIL_REAL result;

  logminusx_re = 0.5L * TSIL_LOG (xsquared);
  t1 = TSIL_dilog_series_real (xinv);
  t2 = 0.5L * (logminusx_re * logminusx_re);
  result = -t1 - t2;

  if (x > 0.0L)
    result += 2.L*zeta2L;
  else
    result += -zeta2L;

  return result;
}

/* ************************************************************** */
/* 
   Uses the Cohen, Lewin and Zagier formula above, to efficiently
   compute the dilog when |x| is not very small or very large.  The
   actual range of absolute convergence is e^-(Pi Sqrt[3]) < |x| <
   e^(Pi Sqrt[3]) The implementation below must not be used when 
   x < 0.
*/

TSIL_REAL TSIL_dilog_CLZseries_real (TSIL_REAL x)
{
  TSIL_REAL u, usquared, ucubed, ufifth;
  TSIL_REAL logminusu_re;
  TSIL_REAL term1, term2, term3, term4;
  TSIL_REAL first4terms;
  TSIL_REAL utothek;
  TSIL_REAL restofterms = 0.0L;
  TSIL_REAL term;
  TSIL_REAL result;
  int j;

  /* 
   Below is an array of precomputed coefficients, needed for the series:

   Li_2(e^u) = zeta2 + u (1- ln(-u)) -(u^2)/4 - (u^3)/72 + (u^5)/14400
               + \sum_{n >= 7} zeta(2-n) (u^n)/n!

   found in Proposition 1 of Cohen, Lewin, and Zagier, Experimental
   Mathematics, volume 1, (1992), p. 26. The first entry in the array
   below is the n=7 term.
  */

  TSIL_REAL CLZcoeffs[] = {
    -7.87351977828168304358781e-7L,
    1.14822163433274544385655e-8L,
    -1.89788699889709990720092e-10L,
    3.38730137095352127233826e-12L,
    -6.37263644318318039658428e-14L,
    1.24620599129506723045228e-15L,
    -2.51054446089995455091693e-17L,
    5.17825880609062350724171e-19L,
    -1.0887357368300848844274e-20L,
    2.32574411430208722345128e-22L,
    -5.03519521314738956081657e-24L,
    1.10264992943812153330081e-25L,
    -2.43865855090073447345264e-27L,
    5.4401426788562523155908e-29L,
    -1.22283401312173521165232e-30L,
    2.76726346896795058422056e-32L,
    -6.30009059183201394873992e-34L,
    1.44208683884184752107295e-35L,
    -3.31709399915954280435211e-37L,
    7.66391355792065788742835e-39L,
    -1.77787147338306578734017e-40L,
    4.13960589823413734492671e-42L,
    -9.67155703608110179257412e-44L,
    2.26671870167661237051842e-45L,
    -5.32795631132825397222586e-47L
  };

  u = 0.5L*TSIL_LOG (x*x);
  usquared = u*u;
  ucubed = usquared*u;
  ufifth = ucubed*usquared;
  logminusu_re = 0.5L * TSIL_LOG (usquared);

  term1 = zeta2L + u - u * logminusu_re;
  term2 = -0.25L * usquared;
  term3 = -ucubed / 72.L;
  term4 = ufifth/14400.L;

  first4terms = term4 + term3 + term2 + term1; 
/*    first4terms += term3; */
/*    first4terms += term2; */
/*    first4terms += term1; */

  utothek = ufifth;

  /* In the following loop, k = 2j+7 */
  for (j = 0; j < 24; j++)
    {
      utothek = utothek * usquared;
      term = utothek*CLZcoeffs[j];
      restofterms += term;
      result = restofterms + first4terms;
      if (TSIL_FABS (term) < 0.1L*TSIL_TOL* TSIL_FABS (result))
	break;
    }

  return result;
}

/* ************************************************************** */
/* 
   We use TSIL_dilog_neg_real to avoid having to use the CLZ series for 
   x < 0, where it converges more slowly. It implements the identity:
   Li2(z) = Li2 (z/(z-1)) + ...
*/

TSIL_REAL TSIL_dilog_neg_real (TSIL_REAL x)
{
  TSIL_REAL t1, t2, t3, t4, t5;
  TSIL_REAL log1minusx_re;

  t1 = 1.L - x;
  t2 = t1 * t1;
  log1minusx_re = 0.5L * TSIL_LOG (t2);
  t3 = -0.5L * log1minusx_re * log1minusx_re;

  if (x > 1.0L) 
    t3 += 3.L*zeta2L;

  t4 = x / (x - 1.L);
  t5 = TSIL_dilog_CLZseries_real (t4);

  return t3 - t5;
}

/* ************************************************************** */

TSIL_COMPLEX TSIL_Dilog (TSIL_COMPLEX z)
{
  TSIL_REAL z_re;
  TSIL_REAL z_im;
  TSIL_REAL res_re, res_im;
  TSIL_REAL absz2; 

  z_re = TSIL_CREAL(z);
  z_im = TSIL_CIMAG(z);
  absz2 = z_re * z_re + z_im * z_im;

  /* First trap the case of real z. */
  if (z_im == 0.0L)
    {
      res_re = TSIL_dilog_real (z_re);
      if (z_re > 1.0L)
        res_im = -M_PI_LONG * TSIL_LOG (z_re);
      else
        res_im = 0.0L;
    }
  else if (absz2 < dilog_crossover_low)
    TSIL_dilog_series_complex (z_re, z_im, &res_re, &res_im);
  else if (absz2 > dilog_crossover_high)
    TSIL_dilog_largeabs_complex (z_re, z_im, &res_re, &res_im);
  else if (z_re < 0)
    TSIL_dilog_neg_complex (z_re, z_im, &res_re, &res_im);
  else
    TSIL_dilog_CLZseries_complex (z_re, z_im, &res_re, &res_im);
  
  return res_re + I*res_im;
}

/* ************************************************************** */
/* Assumes |z| < 1, and will only be called here for |z| < 0.2    */

int 
TSIL_dilog_series_complex (TSIL_REAL z_re,
		      TSIL_REAL z_im, 
		      TSIL_REAL *result_re,
		      TSIL_REAL *result_im)
{
  TSIL_REAL ztothek_re, ztothek_retemp, ztothek_im;
  TSIL_REAL zsquared_re, zsquared_im, zcubed_re, zcubed_im;
  TSIL_REAL term_re, term_im;
  TSIL_REAL ksquared;
  TSIL_REAL absz = TSIL_SQRT (z_re * z_re + z_im * z_im);
  TSIL_REAL res_re = 0.0L;
  TSIL_REAL res_im = 0.0L;
  int kmax, k;

  zsquared_re = z_re * z_re - z_im * z_im;
  zsquared_im = 2.0L * z_re * z_im;
  zcubed_re = zsquared_re * z_re - zsquared_im * z_im;
  zcubed_im = zsquared_re * z_im + zsquared_im * z_re;
  ztothek_re = zcubed_re;
  ztothek_im = zcubed_im;

  kmax = ceil( (float) (TSIL_LOG(TSIL_TOL)/TSIL_LOG(absz)) ) + 1;
  if (kmax < 5)
    kmax=5;

/* 
   Skip the first three terms for now, while adding up the smaller
   terms, to avoid roundoff errors.
*/
  for (k = 4; k < kmax; k++)
    {
      ztothek_retemp = ztothek_re * z_re - ztothek_im * z_im;
      ztothek_im = ztothek_re * z_im + ztothek_im * z_re;
      ztothek_re = ztothek_retemp;
      ksquared = (TSIL_REAL) (k * k);
      term_re = ztothek_re / ksquared;
      term_im = ztothek_im / ksquared;
      res_re += term_re;
      res_im += term_im;
    }

  res_re += zcubed_re / 9.0L;
  res_re += 0.25L * zsquared_re;
  res_re += z_re;
  res_im += zcubed_im / 9.0L;
  res_im += 0.25L * zsquared_im;
  res_im += z_im;

  *result_re = res_re;
  *result_im = res_im;
  return 1;
}

/* ************************************************************** */
/* 
   This reflects z outside the unit circle in the complex plane to
   inside. It implements the identity:
   Li_2 (z) = -Li_2 (1/z) - zeta(2) - (1/2) (log(-z))^2 
*/

int
TSIL_dilog_largeabs_complex (TSIL_REAL z_re,
			TSIL_REAL z_im, 
			TSIL_REAL *result_re,
			TSIL_REAL *result_im)
{
  TSIL_REAL logminusz_re, logminusz_im;
  TSIL_REAL t1_re, t1_im;
  TSIL_REAL zabs2, zin_re, zin_im;

  zabs2 = z_re * z_re + z_im * z_im;
  logminusz_re = 0.5L * TSIL_LOG (zabs2);
  logminusz_im = TSIL_ATAN2 (-z_im, -z_re);
  zin_re = z_re / zabs2;
  zin_im = -z_im / zabs2;

  TSIL_dilog_series_complex (zin_re, zin_im, &t1_re, &t1_im);

  *result_re = -t1_re + 0.5L * (logminusz_im * logminusz_im -
				logminusz_re * logminusz_re) - zeta2L;

  *result_im = -t1_im - logminusz_re * logminusz_im;

  return 1;
}

/* ************************************************************** */

int
TSIL_dilog_neg_complex (TSIL_REAL z_re,
		   TSIL_REAL z_im,
		   TSIL_REAL *result_re,
		   TSIL_REAL *result_im)
{
  TSIL_REAL t1, t2, t3_re, t3_im, t4_re, t4_im, t5_re, t5_im;
  TSIL_REAL log1minusz_re, log1minusz_im;

  t1 = 1.L - z_re;
  t2 = t1 * t1 + z_im * z_im;
  log1minusz_re = 0.5L * TSIL_LOG (t2);
  log1minusz_im = TSIL_ATAN2 (-z_im, t1);
  t3_re = 0.5L * (log1minusz_im * log1minusz_im - log1minusz_re * log1minusz_re);
  t3_im = -log1minusz_re * log1minusz_im;
  t4_re = 1.L + (z_re - 1.L) / t2;
  t4_im = -z_im / t2;
  TSIL_dilog_CLZseries_complex (t4_re, t4_im, &t5_re, &t5_im);

  *result_re = t3_re - t5_re;
  *result_im = t3_im - t5_im;

  return 1;
}

/* ************************************************************** */
/* See earlier comment for TSIL_dilog_CLZseries_real. */

int
TSIL_dilog_CLZseries_complex (TSIL_REAL z_re,
			 TSIL_REAL z_im, 
			 TSIL_REAL *result_re,
			 TSIL_REAL *result_im)
{
  TSIL_REAL u_re, u_im, u_resq, u_imsq;
  TSIL_REAL usquared_re, usquared_im;
  TSIL_REAL ucubed_re, ucubed_im;
  TSIL_REAL ufifth_re, ufifth_im;
  TSIL_REAL logminusu_re, logminusu_im;
  TSIL_REAL term1_re, term2_re, term3_re, term4_re;
  TSIL_REAL term1_im, term2_im, term3_im, term4_im;
  TSIL_REAL first4terms_re, first4terms_im;
  TSIL_REAL utothek_re, utothek_retemp, utothek_im;
  TSIL_REAL restofterms_re = 0.0L;
  TSIL_REAL restofterms_im = 0.0L;
  TSIL_REAL term_re, term_im;
  TSIL_REAL coeffk;
  TSIL_REAL res_re, res_im;
  int j;

  TSIL_REAL CLZcoeffs[] = {
    -7.87351977828168304358781e-7L,
    1.14822163433274544385655e-8L,
    -1.89788699889709990720092e-10L,
    3.38730137095352127233826e-12L,
    -6.37263644318318039658428e-14L,
    1.24620599129506723045228e-15L,
    -2.51054446089995455091693e-17L,
    5.17825880609062350724171e-19L,
    -1.0887357368300848844274e-20L,
    2.32574411430208722345128e-22L,
    -5.03519521314738956081657e-24L,
    1.10264992943812153330081e-25L,
    -2.43865855090073447345264e-27L,
    5.4401426788562523155908e-29L,
    -1.22283401312173521165232e-30L,
    2.76726346896795058422056e-32L,
    -6.30009059183201394873992e-34L,
    1.44208683884184752107295e-35L,
    -3.31709399915954280435211e-37L,
    7.66391355792065788742835e-39L,
    -1.77787147338306578734017e-40L,
    4.13960589823413734492671e-42L,
    -9.67155703608110179257412e-44L,
    2.26671870167661237051842e-45L,
    -5.32795631132825397222586e-47L
  };

  u_re = 0.5L * TSIL_LOG (z_re * z_re + z_im * z_im);
  u_im = TSIL_ATAN2 (z_im, z_re);
  u_resq = u_re * u_re;
  u_imsq = u_im * u_im;

  logminusu_re = 0.5L * TSIL_LOG (u_resq + u_imsq);
  logminusu_im = TSIL_ATAN2 (-u_im, -u_re);

  usquared_re = u_resq - u_imsq;
  usquared_im = 2.0L * u_re * u_im;

  ucubed_re = usquared_re * u_re - usquared_im * u_im;
  ucubed_im = usquared_re * u_im + usquared_im * u_re;

  ufifth_re = ucubed_re * usquared_re - ucubed_im * usquared_im;
  ufifth_im = ucubed_re * usquared_im + ucubed_im * usquared_re;

  term1_re = zeta2L + u_re - u_re * logminusu_re + u_im * logminusu_im;
  term1_im = u_im - u_im * logminusu_re - u_re * logminusu_im;

  term2_re = -0.25L * usquared_re;
  term2_im = -0.25L * usquared_im;

  term3_re = -ucubed_re / 72.L;
  term3_im = -ucubed_im / 72.L;

  term4_re = ufifth_re / 14400.L;
  term4_im = ufifth_im / 14400.L;

  first4terms_re = term4_re + term3_re + term2_re + term1_re;
/*    first4terms_re += term3_re; */
/*    first4terms_re += term2_re; */
/*    first4terms_re += term1_re; */

  first4terms_im = term4_im + term3_im + term2_im + term1_im;
/*    first4terms_im += term3_im; */
/*    first4terms_im += term2_im; */
/*    first4terms_im += term1_im; */

  utothek_re = ufifth_re;
  utothek_im = ufifth_im;

  /* In the following loop, k = 2j+7 */
  for (j = 0; j < 24; j++)
    {
      utothek_retemp = utothek_re * usquared_re - utothek_im * usquared_im;
      utothek_im = utothek_re * usquared_im + utothek_im * usquared_re;
      utothek_re = utothek_retemp;
      coeffk = CLZcoeffs[j];
      term_re = coeffk * utothek_re;
      term_im = coeffk * utothek_im;
      restofterms_re += term_re;
      restofterms_im += term_im;
      res_re = restofterms_re + first4terms_re;
      res_im = restofterms_im + first4terms_im;
      if (
	  (TSIL_FABS (term_re) < 0.1L*TSIL_TOL * TSIL_FABS (res_re)) &&
	  (TSIL_FABS (term_im) < 0.1L*TSIL_TOL * TSIL_FABS (res_im)) )
        break;
    }

  *result_re = res_re;
  *result_im = res_im;

  return 1;
}
