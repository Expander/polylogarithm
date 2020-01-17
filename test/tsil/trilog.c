#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <float.h>
#include "tsil.h"

/* Local routines: */
TSIL_COMPLEX TSIL_Trilogoutofunitdisk (TSIL_COMPLEX);
TSIL_COMPLEX TSIL_Trilogunitdisk (TSIL_COMPLEX);
TSIL_COMPLEX TSIL_TrilogCLZseries (TSIL_COMPLEX);
TSIL_COMPLEX TSIL_Trilogseries (TSIL_COMPLEX);
TSIL_COMPLEX TSIL_TrilogregionA (TSIL_COMPLEX);
TSIL_COMPLEX TSIL_TrilogregionB (TSIL_COMPLEX);

#define cZeta2 1.644934066848226436472415166646025189219L
#define cZeta3 1.202056903159594285399738161511449990765L
#define PI_longdouble 3.141592653589793238462643383279502884197L
#define TSIL_Infinity ((1.0L+ 1.0L*I)/0.0L)

#define trilog_CLZseries_radius 1.2L
#define trilog_powerseries_radius 0.786152L

#define TSIL_Warn(domain, msg)

/*
  Descriptions of the regions inside the unit disk:
   "CLZseries" is for |ln(z)| < trilog_CLZseries_radius and |z|<1;
    Otherwise,
   "series" is for |z| <= trilog_powerseries_radius
   "A" is for Re(z) <= 0 and trilog_powerseries_radius < |z| <= 1.
   "B" is for Re(z) > 0 and trilog_powerseries_radius < |z| <= 1.
*/


/* ************************************************************** */

TSIL_COMPLEX TSIL_Trilog (TSIL_COMPLEX z)
{
  TSIL_COMPLEX result;
  
  if (TSIL_CABS(z) > 1.0L)
    result = TSIL_Trilogoutofunitdisk (z);
  else
    result = TSIL_Trilogunitdisk (z);

  return result;
}

/* ************************************************************** */

TSIL_COMPLEX TSIL_Trilogoutofunitdisk (TSIL_COMPLEX z)
{
  TSIL_COMPLEX result;
  TSIL_COMPLEX logminusz = TSIL_CLOG(-z);

  if (TSIL_CREAL(z) > 1.0L && TSIL_CIMAG((complex double) z) == 0.0)
    logminusz = I * PI_longdouble + TSIL_CLOG(TSIL_CREAL (z));

  result = TSIL_Trilogunitdisk (1.0L/z) -
    logminusz*(cZeta2 + logminusz*logminusz/6.0L);

  return result;
}

/* ************************************************************** */

TSIL_COMPLEX TSIL_Trilogunitdisk (TSIL_COMPLEX z)
{
  TSIL_COMPLEX result;
  TSIL_REAL rez = TSIL_CREAL (z);
  TSIL_REAL absz = TSIL_CABS (z);
  TSIL_REAL absimz = TSIL_FABS (TSIL_CIMAG (z));

  if (TSIL_CABS(z - 1.0L) < 2.0L * TSIL_TOL)
    result = cZeta3;
  else if (TSIL_CABS(z) < 2.0L * TSIL_TOL)
    result = 0.0L;
  else if (TSIL_CABS(TSIL_CLOG(z)) < trilog_CLZseries_radius)
    result = TSIL_TrilogCLZseries (z);
  else if (absz <= trilog_powerseries_radius)
    result = TSIL_Trilogseries (z);
  else if (rez <= 0.0L)
    result = TSIL_TrilogregionA (z);
  else if (rez <= absimz)
    result = TSIL_TrilogregionB (z);
  else {
    TSIL_Warn("TSIL_Trilogunitdisk", "trilog function yielding undefined result.");
    result = TSIL_Infinity;
  }

  return result;
}

/* ************************************************************** */
/* 
   The following is based on the formula by Cohen, Lewin and Zagier.
   It should converge for |ln(z)| < 2 Pi, but is used here only for
   |ln(z)| < trilog_CLZseries_radius. It seems to have a branch cut
   problem for real z > 1, but here it is only used within the unit
   disk anyway. */

TSIL_COMPLEX TSIL_TrilogCLZseries(TSIL_COMPLEX z)
{
  TSIL_COMPLEX logz, logzsquared, logztothek, term;
  TSIL_COMPLEX first6terms, remainingterms, result;
  TSIL_REAL accuracygoal;
  int j;

  TSIL_COMPLEX CLZcoeffs_trilog[] = {
    -9.841899722852103804484756865709246661628e-8L,
    1.14822163433274544385655496766607877719e-9L,
    -1.581572499080916589334097751606169114587e-11L,
    2.419500979252515194527327015649980163421e-13L,
    -3.982897776989487747865172909264620022997e-15L,
    6.923366618305929058068209540950658700107e-17L,
    -1.255272230449977275458465709126553674555e-18L,
    2.353754002768465230564411714140603787718e-20L,
    -4.536398903458687018447507089017008298634e-22L,
    8.94516967039264316712031170773304472046e-24L,
    -1.798284004695496271720202471410154260645e-25L,
    3.675499764793738444336047339126740989995e-27L,
    -7.620807971564795229539485009637654782326e-29L,
    1.600041964369485975173763922573256015938e-30L,
    -3.396761147560375587923120605208518524092e-32L,
    7.282272286757764695317256361444326644631e-34L,
    -1.575022647958003487184978939403782607194e-35L,
    3.433540092480589335887972120165270383376e-37L,
    -7.53884999808987000989116127657764248773e-39L,
    1.666068164765360410310510689355203042974e-40L,
    -3.703898902881387056958685332822841606828e-42L,
    8.279211796468274689853419455144944031424e-44L,
    -1.859914814630981113956561297352068541421e-45L,
    4.19762722532705994540447561866591994072e-47L,
    -9.514207698800453521831900594208537399763e-49L
  };

  logz = TSIL_CLOG(z);
  logzsquared = logz*logz;
  logztothek = logzsquared*logzsquared*logzsquared;

  first6terms = logztothek/86400.0L;
  first6terms += -logzsquared*logzsquared/288.0L;
  first6terms += -logzsquared*logz/12.0L;
  first6terms += (0.75L - 0.5L*TSIL_CLOG(-logz))*logzsquared;
  first6terms += cZeta3 + cZeta2*logz; 
  
  accuracygoal = TSIL_TOL* TSIL_CABS (first6terms);
  remainingterms = 0.0L;

  for (j=0; j< 25; j++)
    {
      logztothek = logztothek*logzsquared;
      term = CLZcoeffs_trilog[j]*logztothek;
      remainingterms += term;
      if (TSIL_CABS (term) < accuracygoal) 
	break;
      if (j == 24) {
        TSIL_Warn("TSIL_TrilogCLZseries", "trilog CLZ series converging too slowly.");
      }
    }

  result = remainingterms + first6terms;

  return result;
}

/* ************************************************************** */

TSIL_COMPLEX TSIL_Trilogseries (TSIL_COMPLEX z)
{
  TSIL_REAL absz = TSIL_CABS (z);
  TSIL_REAL logepsilon = TSIL_LOG (TSIL_TOL);
  TSIL_REAL mlogabsz;
  TSIL_COMPLEX sum = z;
  TSIL_COMPLEX ztothek;
  TSIL_COMPLEX term;
  TSIL_COMPLEX kcubed;
  int k, kmax;

  mlogabsz = -TSIL_CLOG (absz);

/*
  The following kmax is hopefully designed to give accuracy to within
  e^logepsilon, with some safety margin built in. Not completely
  tested, but it seems good enough for government work anyway.
*/

  kmax = 5 + (int) (( 6.0 -logepsilon -3.0 * log(-logepsilon)
		      + 3.0 * log (mlogabsz)) / mlogabsz);

  for (k = kmax; k > 1; k--)
    {
      ztothek = TSIL_CPOW (z, k);
      kcubed = k*k*k;
      term = ztothek/kcubed;
      sum += term;
    }

  return sum;
}

/* ************************************************************** */

TSIL_COMPLEX TSIL_TrilogregionA (TSIL_COMPLEX z)
{
  TSIL_COMPLEX result;
  TSIL_COMPLEX log1minusz = TSIL_CLOG (1.0L - z);
  TSIL_COMPLEX logminusz = TSIL_CLOG (-z);

  result = -TSIL_Trilogunitdisk(1.0L/(1.0L - z)) - TSIL_Trilogunitdisk (z/(z - 1.0L)) 
           + log1minusz*(log1minusz*log1minusz/3.0L -
	     0.5L*logminusz*log1minusz - cZeta2) + cZeta3;

  return result;
}

/* ************************************************************** */

TSIL_COMPLEX TSIL_TrilogregionB (TSIL_COMPLEX z)
{
  TSIL_COMPLEX result;

  result = 0.25L * TSIL_Trilogunitdisk (z * z) - TSIL_Trilogunitdisk (-z);

  return result;
}
