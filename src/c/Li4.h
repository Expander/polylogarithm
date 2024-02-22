/* ====================================================================
 * This file is part of Polylogarithm.
 *
 * Polylogarithm is licenced under the MIT License.
 * ==================================================================== */

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/** real polylogarithm with n=4 (tetralogarithm) with double precision */
double li4(double);

/** complex polylogarithm with n=4 (tetralogarithm) with double precision */
double _Complex cli4(double _Complex);

/** complex polylogarithm with n=4 (tetralogarithm) with long double precision */
long double _Complex cli4l(long double _Complex);

#ifdef __cplusplus
}
#endif
