#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/** real polylogarithm with n=2 (dilogarithm) */
double li2(double x);

/** real polylogarithm with n=2 (dilogarithm) with long double precision */
long double li2l(long double x);

/** complex polylogarithm with n=2 (dilogarithm) */
void cli2_(double re, double im, double* res_re, double* res_im);

/** complex polylogarithm with n=2 (dilogarithm) with long double precision */
void cli2l_(long double re, long double im, long double* res_re, long double* res_im);

/** complex polylogarithm with n=3 (trilogarithm) */
void cli3_(double re, double im, double* res_re, double* res_im);

/** complex polylogarithm with n=3 (trilogarithm) with long double precision */
void cli3l_(long double re, long double im, long double* res_re, long double* res_im);

/** complex polylogarithm with n=4 */
void cli4_(double re, double im, double* res_re, double* res_im);

/** complex polylogarithm with n=4 with long double precision */
void cli4l_(long double re, long double im, long double* res_re, long double* res_im);

/** complex polylogarithm with n=5 */
void cli5_(double re, double im, double* res_re, double* res_im);

/** complex polylogarithm with n=5 with long double precision */
void cli5l_(long double re, long double im, long double* res_re, long double* res_im);

/** complex polylogarithm with n=6 */
void cli6_(double re, double im, double* res_re, double* res_im);

/** complex polylogarithm with n=6 with long double precision */
void cli6l_(long double re, long double im, long double* res_re, long double* res_im);

#ifdef ENABLE_FORTRAN

/** real polylogarithm with n=2 (dilogarithm), Fortran implementation */
double li2_fortran(double x);

/** complex polylogarithm with n=2 (dilogarithm), Fortran implementation */
void cli2_fortran(double re, double im, double* res_re, double* res_im);

/** complex polylogarithm with n=3 (trilogarithm), Fortran implementation */
void cli3_fortran(double re, double im, double* res_re, double* res_im);

/** complex polylogarithm with n=4, Fortran implementation */
void cli4_fortran(double re, double im, double* res_re, double* res_im);

#endif

double algorithm_327(double x);

double algorithm_490(double x);

double babar_dilog(double x);

double cephes_dilog(double x);

double cephes_dilog_2(double x);

double hassani_dilog(double x);

double koelbig_dilog(double x);

long double koelbig_dilogl(long double x);

double morris_dilog(double x);

void hollik_dilog(double re, double im, double* res_re, double* res_im);

void sherpa_dilog(double re, double im, double* res_re, double* res_im);

void spheno_dilog(double re, double im, double* res_re, double* res_im);

long double tsil_dilog_real(long double x);

void tsil_dilog_complex(long double re, long double im, long double* res_re, long double* res_im);

void tsil_trilog_complex(long double re, long double im, long double* res_re, long double* res_im);

#ifdef __cplusplus
}
#endif
