#define TSIL_VERSION "1.44"

/* General header file for the user API and all TSIL types.  This file
   must be included in all applications of TSIL.  */

#ifndef _TSIL_H_
#define _TSIL_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Select basic floating point data size: */

#if defined(TSIL_SIZE_DOUBLE)

typedef double          TSIL_REAL;
typedef double _Complex TSIL_COMPLEX;
#define TSIL_EXP        exp
#define TSIL_CEXP       cexp
#define TSIL_LOG        log
#define TSIL_CLOG       clog
#define TSIL_FABS       fabs
#define TSIL_CABS       cabs
#define TSIL_SQRT       sqrt
#define TSIL_CSQRT      csqrt
#define TSIL_POW        pow
#define TSIL_CPOW       cpow
#define TSIL_ATAN       atan
#define TSIL_ATAN2      atan2
#define TSIL_CREAL      creal
#define TSIL_CIMAG      cimag
#define TSIL_CONJ       conj
#define TSIL_EPSILON    DBL_EPSILON
#define TSIL_TOL        100.0*DBL_EPSILON

#else   /* Assume LONG if we get here: */

typedef long double          TSIL_REAL;
typedef long double _Complex TSIL_COMPLEX;
#define TSIL_EXP             expl
#define TSIL_CEXP            cexpl
#define TSIL_LOG             logl
#define TSIL_CLOG            clogl
#define TSIL_FABS            fabsl
#define TSIL_CABS            cabsl
#define TSIL_SQRT            sqrtl
#define TSIL_CSQRT           csqrtl
#define TSIL_POW             powl
#define TSIL_CPOW            cpowl
#define TSIL_ATAN            atanl
#define TSIL_ATAN2           atan2l
#define TSIL_CREAL           creall
#define TSIL_CIMAG           cimagl
#define TSIL_CONJ            conjl
#define TSIL_EPSILON         LDBL_EPSILON
#define TSIL_TOL             1000.0L*LDBL_EPSILON

#endif /* End selection of basic floating point data size. */


/* ========= Datatype Definitions ========= */

/* ======================================== */
/*            B-type function               */
/* ======================================== */
struct TSIL_Btype {

  int          which;
  TSIL_REAL    arg[2];
  TSIL_COMPLEX value;
  TSIL_COMPLEX deriv;

  /* Evolution coefficients: */
  TSIL_REAL B_den[2]; /* Contains Th and Ps */
  TSIL_REAL B_cB[2];  /* Contains Th/2 and Ps/2 */
  TSIL_REAL B_c[2];   /* Contains ABTh and ABps */
};

typedef struct TSIL_Btype TSIL_BTYPE;


/* ======================================== */
/*            S-type function               */
/* ======================================== */
struct TSIL_Stype {

  int          which;
  TSIL_REAL    arg[3];
  TSIL_COMPLEX value;
  TSIL_COMPLEX deriv;
  TSIL_COMPLEX bold[3]; /* Coefficients of 1/eps^n where n=0,1,2 */

  /* Evolution factor: */
  TSIL_REAL S_c;

  /* Pointers to needed data: */
  TSIL_COMPLEX *tval[3];
};

typedef struct TSIL_Stype TSIL_STYPE;


/* ======================================== */
/*            T-type function               */
/* ======================================== */
struct TSIL_Ttype {

  int          which;
  TSIL_REAL    arg[3];
  TSIL_COMPLEX value;
  TSIL_COMPLEX deriv;
  TSIL_COMPLEX bold[3]; /* Coefficients of 1/eps^n where n=0,1,2 */

  /* Evolution coefficients: */
  TSIL_REAL T_den[4]; /* Contains Th, Px, Py, Pz for use in denominators */
  TSIL_REAL cTS_num[4];
  TSIL_REAL cTT1_num[4];
  TSIL_REAL cTT2_num[4];
  TSIL_REAL cTT3_num[4];
  TSIL_REAL cTs_num, cT_num[4]; /* Contains bT, aT, axT, ayT, azT */

  /* Pointers to needed data: */
  TSIL_COMPLEX *sval;
  TSIL_COMPLEX *tval[3];
};

typedef struct TSIL_Ttype TSIL_TTYPE;


/* ======================================== */
/*            U-type function               */
/* ======================================== */
struct TSIL_Utype {

  int          which;
  TSIL_REAL    arg[4];
  TSIL_COMPLEX value;
  TSIL_COMPLEX deriv;
  TSIL_COMPLEX bold[3]; /* Coefficients of 1/eps^n where n=0,1,2 */

  /* Evolution coefficients: */
  TSIL_REAL den_th;    
  TSIL_REAL den_ps;    
  TSIL_REAL cU_th; 
  TSIL_REAL cU_ps; 
  TSIL_REAL cS_th; 
  TSIL_REAL cS_ps; 
  TSIL_REAL cT1_th;
  TSIL_REAL cT1_ps;
  TSIL_REAL cT2;    
  TSIL_REAL cT3;    
  TSIL_REAL con;      

  /* Pointers to needed data: */
  TSIL_COMPLEX *sval;
  TSIL_COMPLEX *tval[3];
};

typedef struct TSIL_Utype TSIL_UTYPE;


/* ======================================== */
/*            V-type function               */
/* ======================================== */
struct TSIL_Vtype {

  int          which;
  int          isSet;
  TSIL_REAL    arg[4];
  TSIL_COMPLEX value;
  TSIL_COMPLEX bold[3]; /* Coefficients of 1/eps^n where n=0,1,2 */
};

typedef struct TSIL_Vtype TSIL_VTYPE;


/* ======================================== */
/*            Tbar-type function            */
/* ======================================== */
struct TSIL_Tbartype {

  int          which;
  int          isSet;
  TSIL_COMPLEX value;
  TSIL_REAL    arg[3];
};

typedef struct TSIL_Tbartype TSIL_TbarTYPE;


/* ======================================== */
/*            M-type function               */
/* ======================================== */
struct TSIL_Mtype {

  TSIL_REAL    arg[5];
  TSIL_COMPLEX value;
  TSIL_COMPLEX deriv;

  /* Evolution coefficients: */
  int       extramassdim1, extramassdim2;
  TSIL_REAL THxz, THyu, PSxz, PSyu;
  TSIL_REAL adenom[3];
  TSIL_REAL aMU[4], aMUs[4], bMU[4][2];
  TSIL_REAL aMT[5], bMT[5][4];
  TSIL_REAL aMT5s;
  TSIL_REAL aMS, bMS[4], cMSconst;
  TSIL_REAL dMB, dMBs;
  TSIL_REAL aMB[2], bMB[2][2];
  TSIL_REAL aM, aMs, bM[4];

  /* Pointers to needed data: */
  TSIL_COMPLEX *bval[2];
  TSIL_COMPLEX *sval[2];
  TSIL_COMPLEX *tval[6];
  TSIL_COMPLEX *uval[4];
};

typedef struct TSIL_Mtype TSIL_MTYPE;

/* ======================================== */
/*          General Data Structure          */
/* ======================================== */
struct TSIL_Data
{
  int isAligned;     /* Indicates whether the data sub-objects have
                        been constructed */
  int isInitialized; /* Indicates whether initial values have been set */

  int status;        /* Current status: evaluated and, if so, how (see
                        tsil_global.h for value definitions) */
  int whichFns;      /* Indicates which functions are (or have) been 
			evaluated (values in tsil_global.h) */

  /* Data */
  TSIL_REAL     x, y, z, u, v, s, qq;
  TSIL_BTYPE    B[2];
  TSIL_STYPE    S[2];
  TSIL_TTYPE    T[6];
  TSIL_TbarTYPE Tbar[6];
  TSIL_UTYPE    U[4];
  TSIL_VTYPE    V[4];
  TSIL_MTYPE    M;

  TSIL_REAL scaleFac;
  TSIL_REAL threshold[4], pseudoThreshold[8];
  TSIL_REAL threshMin, psThreshMin;
  int nThresh, nPthresh;

  /* Initialization point for numerical integration; may be zero or
     SINIT (defined in tsil_params.h) */
  TSIL_COMPLEX sInitial;

  /* Parameters for automatic step-size control in Runge-Kutta */
  TSIL_REAL precisionGoal;
  int       nStepsStart, nStepsMaxCon, nStepsMaxVar, nStepsMin;

  /* Pointers to RK stepper functions: */
  int  (*RKstepper6) ();
  void (*RKstepper5) ();
};

typedef struct TSIL_Data TSIL_DATA;


/* ======================================== */
/*          Results Data Structure          */
/* ======================================== */
struct TSIL_Result {

  TSIL_REAL    x, y, z, u, v, s, qq;
  TSIL_COMPLEX M;
  TSIL_COMPLEX U[4];
  TSIL_COMPLEX V[4];
  TSIL_COMPLEX T[6];
  TSIL_COMPLEX S[2];
  TSIL_COMPLEX B[2];
  TSIL_COMPLEX TBAR[6];
};

typedef struct TSIL_Result TSIL_RESULT;

/* For permutations of TSIL_RESULTs (see function TSIL_PermuteResult) */
enum {NOSWAP, XYandZU, XZandYU, XUandYZ};

/* Toggle to control printing of warning messages */
extern int printTSILWarns;
/* extern FILE *warnfile, *errfile; */

/* === Prototypes for functions in the user API: === */

/* Basic evaluation functions: */
int TSIL_SetParameters (TSIL_DATA *, TSIL_REAL, TSIL_REAL, TSIL_REAL,
			TSIL_REAL, TSIL_REAL, TSIL_REAL);
int TSIL_SetParametersSTU (TSIL_DATA *, TSIL_REAL, TSIL_REAL,
			   TSIL_REAL, TSIL_REAL, TSIL_REAL);
int TSIL_SetParametersST (TSIL_DATA *, TSIL_REAL, TSIL_REAL, 
			  TSIL_REAL, TSIL_REAL);
int TSIL_Evaluate (TSIL_DATA *, TSIL_REAL);
int TSIL_GetStatus (TSIL_DATA *);
int TSIL_GetData (TSIL_DATA *, const char *, TSIL_COMPLEX *);
int TSIL_GetDataR (TSIL_RESULT *, const char *, TSIL_COMPLEX *);
int TSIL_GetBoldData (TSIL_DATA *, const char *, TSIL_COMPLEX [][3]);
TSIL_COMPLEX TSIL_GetFunction (TSIL_DATA *, const char *);
TSIL_COMPLEX TSIL_GetFunctionR (TSIL_RESULT *, const char *);
TSIL_COMPLEX TSIL_GetBoldFunction (TSIL_DATA *, const char *, int);
int TSIL_NumFuncs (const char *);

/* I/O and related functions: */
void TSIL_PrintStatus (TSIL_DATA *);
void TSIL_PrintData   (TSIL_DATA *);
void TSIL_WriteData   (FILE *, TSIL_DATA *);
void TSIL_PrintDataM  (TSIL_DATA *);
void TSIL_WriteDataM  (FILE *, TSIL_DATA *);
void TSIL_cprintf     (TSIL_COMPLEX);
void TSIL_cprintfM    (TSIL_COMPLEX);
void TSIL_Error       (const char *, const char *, int);
void TSIL_Warn        (const char *, const char *);
void TSIL_WarnsOn     (void);
void TSIL_WarnsOff    (void);
void TSIL_PrintInfo   (void);
void TSIL_Info        (const char *);

/* Utilities: */
int  TSIL_CopyResult (TSIL_DATA *, TSIL_RESULT *);
void TSIL_PermuteResult (TSIL_RESULT *, int, TSIL_RESULT *);
void TSIL_ResetStepSizeParams (TSIL_DATA *,TSIL_REAL, int, int, int, int);
int  TSIL_IsInfinite (TSIL_COMPLEX);
int  TSIL_DataSize (void);
void TSIL_PrintVersion (void);

/* Analytic cases: */
TSIL_COMPLEX TSIL_Dilog  (TSIL_COMPLEX);
TSIL_COMPLEX TSIL_Trilog (TSIL_COMPLEX);
TSIL_COMPLEX TSIL_A      (TSIL_REAL, TSIL_REAL);
TSIL_COMPLEX TSIL_Ap     (TSIL_REAL, TSIL_REAL);
TSIL_COMPLEX TSIL_Aeps   (TSIL_REAL, TSIL_REAL);
TSIL_COMPLEX TSIL_B      (TSIL_REAL, TSIL_REAL, TSIL_COMPLEX, TSIL_REAL);
TSIL_COMPLEX TSIL_Bp     (TSIL_REAL, TSIL_REAL, TSIL_COMPLEX, TSIL_REAL);
TSIL_COMPLEX TSIL_dBds   (TSIL_REAL, TSIL_REAL, TSIL_COMPLEX, TSIL_REAL);
TSIL_COMPLEX TSIL_Beps   (TSIL_REAL, TSIL_REAL, TSIL_COMPLEX, TSIL_REAL);
TSIL_COMPLEX TSIL_I2     (TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL);
TSIL_COMPLEX TSIL_I2p    (TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL);
TSIL_COMPLEX TSIL_I2p2   (TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL);
TSIL_COMPLEX TSIL_I2pp   (TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL);
TSIL_COMPLEX TSIL_I2p3   (TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL);
int TSIL_Sanalytic    (TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_COMPLEX,
		       TSIL_REAL, TSIL_COMPLEX *);
int TSIL_Tanalytic    (TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_COMPLEX,
		       TSIL_REAL, TSIL_COMPLEX *);
int TSIL_Tbaranalytic (TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_COMPLEX,
		       TSIL_REAL, TSIL_COMPLEX *);
int TSIL_Uanalytic    (TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL,
		       TSIL_COMPLEX, TSIL_REAL, TSIL_COMPLEX *);
int TSIL_Vanalytic    (TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL,
		       TSIL_COMPLEX, TSIL_REAL, TSIL_COMPLEX *);
int TSIL_Manalytic    (TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL,
		       TSIL_REAL, TSIL_COMPLEX, TSIL_COMPLEX *);

#ifdef __cplusplus
}
#endif

#endif /* tsil.h */
