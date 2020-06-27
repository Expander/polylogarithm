/* ====================================================================
 * This file is part of Polylogarithm.
 *
 * Polylogarithm is licenced under the GNU Lesser General Public
 * License (GNU LGPL) version 3.
 * ==================================================================== */

#pragma once

#include <complex.h>
#include <math.h>


static inline double _Complex fast_clog(double _Complex z)
{
   const double rz = creal(z);
   const double iz = cimag(z);
   double arg = atan2(iz, rz);

   if (iz == 0.0 && arg != 0.0)
      arg = -arg;

   return 0.5*log(rz*rz + iz*iz) + I*arg;

   /* const double rz = creal(z) == 0.0 ? fabs(creal(z)) : creal(z); */
   /* const double iz = cimag(z) == 0.0 ? fabs(cimag(z)) : cimag(z); */

   /* return 0.5*log(rz*rz + iz*iz) + I*atan2(iz, rz); */
}


static inline long double _Complex fast_clogl(long double _Complex z)
{
   const long double rz = creall(z) == 0.0L ? fabsl(creall(z)) : creall(z);
   const long double iz = cimagl(z) == 0.0L ? fabsl(cimagl(z)) : cimagl(z);

   return 0.5L*logl(rz*rz + iz*iz) + I*atan2l(iz, rz);
}

