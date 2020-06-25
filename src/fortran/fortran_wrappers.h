#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#ifdef ENABLE_FORTRAN

/** real polylogarithm with n=2 (dilogarithm), Fortran implementation */
void dli2_fortran(const double* x, double* res);

/** complex polylogarithm with n=2 (dilogarithm), Fortran implementation */
void cdli2_fortran(const double* re, const double* im, double* res_re, double* res_im);

/** complex polylogarithm with n=3 (trilogarithm), Fortran implementation */
void cdli3_fortran(const double* re, const double* im, double* res_re, double* res_im);

/** complex polylogarithm with n=4, Fortran implementation */
void cdli4_fortran(const double* re, const double* im, double* res_re, double* res_im);

#endif

#ifdef __cplusplus
}
#endif
