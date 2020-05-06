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

/** complex polylogarithm with n=3 (trilogarithm) */
void cli3_(double re, double im, double* res_re, double* res_im);

#ifdef __cplusplus
}
#endif
