/* ====================================================================
 * This file is part of Polylogarithm.
 *
 * Polylogarithm is licenced under the MIT License.
 * ==================================================================== */

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/** real polylogarithm with n=3 (trilogarithm) with double precision */
double li3(double);

/** complex polylogarithm with n=3 (trilogarithm) with double precision */
double _Complex cli3(double _Complex);

/** complex polylogarithm with n=3 (trilogarithm) with long double precision */
long double _Complex cli3l(long double _Complex);

#ifdef __cplusplus
}
#endif
