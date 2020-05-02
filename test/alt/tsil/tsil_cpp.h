/* TSIL v1.44 */

/* General header file for the user API and all TSIL types.  This file
   should be included in all C++ applications wishing to use TSIL. */

#ifndef TSIL_CPP_H
#define TSIL_CPP_H

#include "tsil.h"
#include <complex>

#ifdef complex
#undef complex
#endif

/* Select basic floating point data size: */

#if defined(TSIL_SIZE_DOUBLE)
typedef std::complex<double> TSIL_COMPLEXCPP;
#ifndef I
#define I            (std::complex<double>(0.0, 1.0))
#endif

#else   /* Assume LONG if we get here: */
typedef std::complex<long double> TSIL_COMPLEXCPP;
#ifndef I
#define I            (std::complex<long double>(0.0L, 1.0L))
#endif

#endif /* End selection of basic floating point data size. */

// ==================================================================
//    C++ wrappers
// ==================================================================
// General wrapper for C functions returning a complex:
inline TSIL_COMPLEXCPP c2cpp (const TSIL_COMPLEX& cval)
{
  return *((TSIL_COMPLEXCPP*) (&cval));
}

// Specific wrappers for existing TSIL functions:
// ==================================================================
inline int TSIL_SetParameters_ (TSIL_DATA *foo, TSIL_REAL x, TSIL_REAL y, TSIL_REAL z,
			 TSIL_REAL u, TSIL_REAL v, TSIL_REAL qq)
{
  return TSIL_SetParameters (foo, x, y, z, u, v, qq);
}
// ==================================================================
inline int TSIL_SetParametersSTU_ (TSIL_DATA *foo, TSIL_REAL x, TSIL_REAL y,
			    TSIL_REAL z, TSIL_REAL u, TSIL_REAL qq)
{
  return TSIL_SetParametersSTU (foo, x, y, z, u, qq);
}
// ==================================================================
inline int TSIL_SetParametersST_ (TSIL_DATA *foo, TSIL_REAL x, TSIL_REAL y, 
			  TSIL_REAL z, TSIL_REAL qq)
{
  return TSIL_SetParametersST (foo, x, y, z, qq);
}
// ==================================================================
inline int TSIL_Evaluate_ (TSIL_DATA * foo, TSIL_REAL ss)
{
  return TSIL_Evaluate (foo, ss);
}
// ==================================================================
inline int TSIL_GetStatus_ (TSIL_DATA *foo)
{
  return TSIL_GetStatus (foo);
}
// ==================================================================
inline int TSIL_GetData_ (TSIL_DATA *foo, const char *which, TSIL_COMPLEXCPP *res)
{
  return TSIL_GetData (foo, which, (TSIL_COMPLEX *) res);
}
// ==================================================================
inline int TSIL_GetBoldData_ (TSIL_DATA *foo, const char *which, TSIL_COMPLEXCPP res[][3])
{
  TSIL_COMPLEX scratch[6][3];
  int i, j, retval, nfuncs = TSIL_NumFuncs (which);

  retval = TSIL_GetBoldData (foo, which, scratch);
  for (i=0; i<nfuncs; i++)
    for (j=0; j<3; j++)
      res[i][j] = c2cpp (scratch[i][j]);

  return retval;
}
// ==================================================================
inline TSIL_COMPLEXCPP TSIL_GetFunction_ (TSIL_DATA *foo, const char *which)
{
  return c2cpp (TSIL_GetFunction (foo, which));
}
// ==================================================================
inline TSIL_COMPLEXCPP TSIL_GetBoldFunction_ (TSIL_DATA *foo, const char *which, int m)
{
  return c2cpp (TSIL_GetBoldFunction (foo, which, m));
}
// ==================================================================
/* I/O and related functions: */
inline void TSIL_PrintStatus_ (TSIL_DATA *foo)
{
  TSIL_PrintStatus (foo);
}
// ==================================================================
inline void TSIL_PrintData_ (TSIL_DATA *foo)
{
  TSIL_PrintData (foo);
}
// ==================================================================
inline void TSIL_WriteData_ (FILE *fp, TSIL_DATA *foo)
{
  TSIL_WriteData (fp, foo);
}
// ==================================================================
inline void TSIL_PrintDataM_ (TSIL_DATA *foo)
{
  TSIL_PrintDataM (foo);
}
// ==================================================================
inline void TSIL_WriteDataM_ (FILE *fp, TSIL_DATA *foo)
{
  TSIL_WriteData (fp, foo);
}
// ==================================================================
inline void TSIL_cprintf_ (TSIL_COMPLEXCPP z)
{
  TSIL_COMPLEX zz = *((TSIL_COMPLEX*) (&z));
  TSIL_cprintf (zz);
}
// ==================================================================
inline void TSIL_cprintfM_ (TSIL_COMPLEXCPP z)
{
  TSIL_COMPLEX zz = *((TSIL_COMPLEX*) (&z));
  TSIL_cprintfM (zz);
}
// ==================================================================
inline void TSIL_Error_ (const char *foo, const char *bar, int m)
{
  TSIL_Error (foo, bar, m);
}
// ==================================================================
inline void TSIL_Warn_ (char *foo, char *bar)
{
  TSIL_Warn (foo, bar);
}
// ==================================================================
inline void TSIL_PrintInfo_ (void)
{
  TSIL_PrintInfo ();
}
// ==================================================================
inline void TSIL_Info_ (char *msg)
{
  TSIL_Info (msg);
}
// ==================================================================
inline int  TSIL_CopyResult_ (TSIL_DATA *foo, TSIL_RESULT *bar)
{
  return TSIL_CopyResult (foo, bar);
}
// ==================================================================
inline void TSIL_PermuteResult_ (TSIL_RESULT *foo, int which, TSIL_RESULT *bar)
{
  TSIL_PermuteResult (foo, which, bar);
}
// ==================================================================
inline void TSIL_ResetStepSizeParams_ (TSIL_DATA *foo,TSIL_REAL z, int i, int j, int k, int l)
{
  TSIL_ResetStepSizeParams (foo, z, i, j, k, l);
}
// ==================================================================
inline int  TSIL_IsInfinite_ (TSIL_COMPLEXCPP z)
{
  TSIL_COMPLEX zz = *((TSIL_COMPLEX*) (&z));
  return TSIL_IsInfinite (zz);
}
// ==================================================================
inline int  TSIL_DataSize_ (void)
{
  return TSIL_DataSize ();
}
// ==================================================================
inline void TSIL_PrintVersion_ (void)
{
  TSIL_PrintVersion ();
}
// ==================================================================
inline TSIL_COMPLEXCPP TSIL_Dilog_ (TSIL_COMPLEXCPP z)
{
  TSIL_COMPLEX zz = *((TSIL_COMPLEX*) (&z));
  return c2cpp (TSIL_Dilog (zz));
}
// ==================================================================
inline TSIL_COMPLEXCPP TSIL_Trilog_ (TSIL_COMPLEXCPP z)
{
  TSIL_COMPLEX zz = *((TSIL_COMPLEX*) (&z));
  return c2cpp (TSIL_Trilog (zz));
}
// ==================================================================
inline TSIL_COMPLEXCPP TSIL_A_ (TSIL_REAL x, TSIL_REAL qq) 
{
  return c2cpp (TSIL_A (x, qq));
}
// ==================================================================
inline TSIL_COMPLEXCPP TSIL_Ap_ (TSIL_REAL x, TSIL_REAL qq)
{
  return c2cpp (TSIL_Ap (x, qq));
}
// ==================================================================
inline TSIL_COMPLEXCPP TSIL_Aeps_ (TSIL_REAL x, TSIL_REAL qq)
{
  return c2cpp (TSIL_Aeps (x, qq));
}
// ==================================================================
inline TSIL_COMPLEXCPP TSIL_B_ (TSIL_REAL x, TSIL_REAL y, TSIL_COMPLEXCPP s, TSIL_REAL qq)
{
  TSIL_COMPLEX ss = *((TSIL_COMPLEX*) (&s));

  return c2cpp (TSIL_B (x, y, ss, qq));
}
// ==================================================================
inline TSIL_COMPLEXCPP TSIL_Bp_ (TSIL_REAL x, TSIL_REAL y, TSIL_COMPLEXCPP s, TSIL_REAL qq)
{
  TSIL_COMPLEX ss = *((TSIL_COMPLEX*) (&s));

  return c2cpp (TSIL_Bp (x, y, ss, qq));
}
// ==================================================================
inline TSIL_COMPLEXCPP TSIL_dBds_ (TSIL_REAL x, TSIL_REAL y, TSIL_COMPLEXCPP s, TSIL_REAL qq)
{
  TSIL_COMPLEX ss = *((TSIL_COMPLEX*) (&s));

  return c2cpp (TSIL_dBds (x, y, ss, qq));
}
// ==================================================================
inline TSIL_COMPLEXCPP TSIL_Beps_ (TSIL_REAL x, TSIL_REAL y, TSIL_COMPLEXCPP s, TSIL_REAL qq)
{
  TSIL_COMPLEX ss = *((TSIL_COMPLEX*) (&s));

  return c2cpp (TSIL_Beps (x, y, ss, qq));
}
// ==================================================================
inline TSIL_COMPLEXCPP TSIL_I2_ (TSIL_REAL x, TSIL_REAL y, TSIL_REAL z, TSIL_REAL qq)
{
  return c2cpp (TSIL_I2 (x, y, z, qq));
}
// ==================================================================
inline TSIL_COMPLEXCPP TSIL_I2p_ (TSIL_REAL x, TSIL_REAL y, TSIL_REAL z, TSIL_REAL qq)
{
  return c2cpp (TSIL_I2p (x, y, z, qq));
}
// ==================================================================
inline TSIL_COMPLEXCPP TSIL_I2p2_ (TSIL_REAL x, TSIL_REAL y, TSIL_REAL z, TSIL_REAL qq)
{
  return c2cpp (TSIL_I2p2 (x, y, z, qq));
}
// ==================================================================
inline TSIL_COMPLEXCPP TSIL_I2pp_ (TSIL_REAL x, TSIL_REAL y, TSIL_REAL z, TSIL_REAL qq)
{
  return c2cpp (TSIL_I2pp (x, y, z, qq));
}
// ==================================================================
inline TSIL_COMPLEXCPP TSIL_I2p3_ (TSIL_REAL x, TSIL_REAL y, TSIL_REAL z, TSIL_REAL qq)
{
  return c2cpp (TSIL_I2p3 (x, y, z, qq));
}
// ==================================================================
// We could just cast the return pointer, e.g.,
//    return TSIL_Sanalytic (x, y, z, ss, qq, (TSIL_COMPLEX *) res);
// (which seems to work) but this seems safer albeit more complicated:

inline int TSIL_Sanalytic_ (TSIL_REAL x, TSIL_REAL y, TSIL_REAL z, TSIL_COMPLEXCPP s,
		     TSIL_REAL qq, TSIL_COMPLEXCPP *res)
{
  TSIL_COMPLEX ss = *((TSIL_COMPLEX*) (&s));
  TSIL_COMPLEX cres;
  int retval;

  retval = TSIL_Sanalytic (x, y, z, ss, qq, &cres);
  *res = *((TSIL_COMPLEXCPP *) (&cres));
  return retval;
}
// ==================================================================
inline int TSIL_Tanalytic_ (TSIL_REAL x, TSIL_REAL y, TSIL_REAL z, TSIL_COMPLEXCPP s,
		     TSIL_REAL qq, TSIL_COMPLEXCPP *res)
{
  TSIL_COMPLEX ss = *((TSIL_COMPLEX*) (&s));
  TSIL_COMPLEX cres;
  int retval;

  retval = TSIL_Tanalytic (x, y, z, ss, qq, &cres);
  *res = *((TSIL_COMPLEXCPP *) (&cres));
  return retval;

  /* return TSIL_Tanalytic (x, y, z, ss, qq, (TSIL_COMPLEX *) res); */
}
// ==================================================================
inline int TSIL_Tbaranalytic_ (TSIL_REAL x, TSIL_REAL y, TSIL_REAL z, TSIL_COMPLEXCPP s,
			TSIL_REAL qq, TSIL_COMPLEXCPP *res)
{
  TSIL_COMPLEX ss = *((TSIL_COMPLEX*) (&s));
  TSIL_COMPLEX cres;
  int retval;

  retval = TSIL_Tbaranalytic (x, y, z, ss, qq, &cres);
  *res = *((TSIL_COMPLEXCPP *) (&cres));
  return retval;

  /* return TSIL_Tbaranalytic (x, y, z, ss, qq, (TSIL_COMPLEX *) res); */
}
// ==================================================================
inline int TSIL_Uanalytic_ (TSIL_REAL x, TSIL_REAL y, TSIL_REAL z, TSIL_REAL u,
		     TSIL_COMPLEXCPP s, TSIL_REAL qq, TSIL_COMPLEXCPP *res)
{
  TSIL_COMPLEX ss = *((TSIL_COMPLEX*) (&s));
  TSIL_COMPLEX cres;
  int retval;

  retval = TSIL_Uanalytic (x, y, z, u, ss, qq, &cres);
  *res = *((TSIL_COMPLEXCPP *) (&cres));
  return retval;

  /* return TSIL_Uanalytic (x, y, z, u, ss, qq, (TSIL_COMPLEX *) res); */
}
// ==================================================================
inline int TSIL_Vanalytic_ (TSIL_REAL x, TSIL_REAL y, TSIL_REAL z, TSIL_REAL u,
		     TSIL_COMPLEXCPP s, TSIL_REAL qq, TSIL_COMPLEXCPP *res)
{
  TSIL_COMPLEX ss = *((TSIL_COMPLEX*) (&s));
  TSIL_COMPLEX cres;
  int retval;

  retval = TSIL_Vanalytic (x, y, z, u, ss, qq, &cres);
  *res = *((TSIL_COMPLEXCPP *) (&cres));
  return retval;

  /* return TSIL_Vanalytic (x, y, z, u, ss, qq, (TSIL_COMPLEX *) res); */
}
// ==================================================================
inline int TSIL_Manalytic_ (TSIL_REAL x, TSIL_REAL y, TSIL_REAL z, TSIL_REAL u, TSIL_REAL v,
		     TSIL_COMPLEXCPP s, TSIL_COMPLEXCPP *res)
{
  TSIL_COMPLEX ss = *((TSIL_COMPLEX*) (&s));
  TSIL_COMPLEX cres;
  int retval;

  retval = TSIL_Manalytic (x, y, z, u, v, ss, &cres);
  *res = *((TSIL_COMPLEXCPP *) (&cres));
  return retval;
}

#endif /* tsil.h */
