/* ====================================================================
 * This file is part of Polylogarithm.
 *
 * Polylogarithm is licenced under the GNU Lesser General Public
 * License (GNU LGPL) version 3.
 * ==================================================================== */

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/** real polylogarithm with n=2 (dilogarithm) */
double li2(double x);

/** complex polylogarithm with n=2 (dilogarithm) */
void cli2_(double re, double im, double* res_re, double* res_im);

#ifdef __cplusplus
}
#endif
