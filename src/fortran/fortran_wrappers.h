#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#ifdef ENABLE_FORTRAN

/** Clausen function with n=2, Fortran implementation */
void cl2_fortran(const double* x, double* res);

/** Clausen function with n=3, Fortran implementation */
void cl3_fortran(const double* x, double* res);

/** Clausen function with n=4, Fortran implementation */
void cl4_fortran(const double* x, double* res);

/** Clausen function with n=5, Fortran implementation */
void cl5_fortran(const double* x, double* res);

/** Clausen function with n=6, Fortran implementation */
void cl6_fortran(const double* x, double* res);

/** real polylogarithm with n=2 (dilogarithm), Fortran implementation */
void li2_fortran(const double* x, double* res);

/** real polylogarithm with n=3 (trilogarithm), Fortran implementation */
void li3_fortran(const double* x, double* res);

/** complex polylogarithm with n=2 (dilogarithm), Fortran implementation */
void cli2_fortran(const double* re, const double* im, double* res_re, double* res_im);

/** complex polylogarithm with n=3 (trilogarithm), Fortran implementation */
void cli3_fortran(const double* re, const double* im, double* res_re, double* res_im);

/** complex polylogarithm with n=4, Fortran implementation */
void cli4_fortran(const double* re, const double* im, double* res_re, double* res_im);

/** complex polylogarithm with n=5, Fortran implementation */
void cli5_fortran(const double* re, const double* im, double* res_re, double* res_im);

/** complex polylogarithm with n=6, Fortran implementation */
void cli6_fortran(const double* re, const double* im, double* res_re, double* res_im);

#endif

#ifdef __cplusplus
}
#endif
