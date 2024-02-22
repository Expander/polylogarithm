/* ====================================================================
 * This file is part of Polylogarithm.
 *
 * Polylogarithm is licenced under the MIT License.
 * ==================================================================== */

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/** real polylogarithm with n=2 (dilogarithm) with single precision */
float li2f(float);

/** real polylogarithm with n=2 (dilogarithm) with double precision */
double li2(double);

/** real polylogarithm with n=2 (dilogarithm) with long double precision */
long double li2l(long double);

/** complex polylogarithm with n=2 (dilogarithm) with single precision */
float _Complex cli2f(float _Complex);

/** complex polylogarithm with n=2 (dilogarithm) with double precision */
double _Complex cli2(double _Complex);

/** complex polylogarithm with n=2 (dilogarithm) with long double precision */
long double _Complex cli2l(long double _Complex);

#ifdef __cplusplus
}
#endif
