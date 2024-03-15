/* ====================================================================
 * This file is part of Polylogarithm.
 *
 * Polylogarithm is licenced under the MIT License.
 * ==================================================================== */

#pragma once

/** returns clog(1 + z) with double precision */
double _Complex clog1p(double _Complex z);

/** returns clog(1 + z) with long double precision */
long double _Complex clog1pl(long double _Complex z);

/** returns clog(z) with double precision */
double _Complex pos_clog(double _Complex z);

/** returns clog(z) with long double precision */
long double _Complex pos_clogl(long double _Complex z);
