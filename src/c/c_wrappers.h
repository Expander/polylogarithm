/* ====================================================================
 * This file is part of Polylogarithm.
 *
 * Polylogarithm is licenced under the MIT License.
 * ==================================================================== */

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/** Clausen function with n=2 */
double cl2(double x);

/** Clausen function with n=3 */
double cl3(double x);

/** Clausen function with n=4 */
double cl4(double x);

/** Clausen function with n=5 */
double cl5(double x);

/** Clausen function with n=6 */
double cl6(double x);

/** Clausen function with n=2 with long double precision */
long double cl2l(long double x);

/** Clausen function with n=3 with long double precision */
long double cl3l(long double x);

/** Clausen function with n=4 with long double precision */
long double cl4l(long double x);

/** Clausen function with n=5 with long double precision */
long double cl5l(long double x);

/** Clausen function with n=6 with long double precision */
long double cl6l(long double x);

/** real polylogarithm with n=2 (dilogarithm) with single precision */
float li2f(float x);

/** real polylogarithm with n=2 (dilogarithm) */
double li2(double x);

/** real polylogarithm with n=3 (trilogarithm) */
double li3(double x);

/** real polylogarithm with n=4 (trilogarithm) */
double li4(double x);

/** real polylogarithm with n=2 (dilogarithm) with long double precision */
long double li2l(long double x);

/** complex polylogarithm with n=2 (dilogarithm) with single precision */
void cli2f_c(float re, float im, float* res_re, float* res_im);

/** complex polylogarithm with n=2 (dilogarithm) with double precision */
void cli2_c(double re, double im, double* res_re, double* res_im);

/** complex polylogarithm with n=2 (dilogarithm) with long double precision */
void cli2l_c(long double re, long double im, long double* res_re, long double* res_im);

/** complex polylogarithm with n=3 (trilogarithm) */
void cli3_c(double re, double im, double* res_re, double* res_im);

/** complex polylogarithm with n=3 (trilogarithm) with long double precision */
void cli3l_c(long double re, long double im, long double* res_re, long double* res_im);

/** complex polylogarithm with n=4 */
void cli4_c(double re, double im, double* res_re, double* res_im);

/** complex polylogarithm with n=4 with long double precision */
void cli4l_c(long double re, long double im, long double* res_re, long double* res_im);

/** complex polylogarithm with n=5 */
void cli5_c(double re, double im, double* res_re, double* res_im);

/** complex polylogarithm with n=5 with long double precision */
void cli5l_c(long double re, long double im, long double* res_re, long double* res_im);

/** complex polylogarithm with n=6 */
void cli6_c(double re, double im, double* res_re, double* res_im);

/** complex polylogarithm with n=6 with long double precision */
void cli6l_c(long double re, long double im, long double* res_re, long double* res_im);

#ifdef __cplusplus
}
#endif
