/* ====================================================================
 * This file is part of Polylogarithm.
 *
 * Polylogarithm is licenced under the MIT License.
 * ==================================================================== */

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/** real polylogarithm with n=2 (dilogarithm) */
double li2(double x);

/** real polylogarithm with n=2 (dilogarithm) with long double precision */
long double li2l(long double x);

/** complex polylogarithm with n=2 (dilogarithm) */
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
