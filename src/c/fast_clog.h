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

   return 0.5*log(rz*rz + iz*iz) + I*atan2(iz, rz);
}


static inline double _Complex fast_pos_clog(double _Complex z)
{
   const double rz = creal(z);
   const double iz = cimag(z);
   double arg = atan2(iz, rz);

   if (iz == 0.0 && arg < 0.0)
      arg = -arg;

   return 0.5*log(rz*rz + iz*iz) + I*arg;
}


static inline long double _Complex fast_clogl(long double _Complex z)
{
   const long double rz = creall(z);
   const long double iz = cimagl(z);

   return 0.5L*logl(rz*rz + iz*iz) + I*atan2l(iz, rz);
}


static inline long double _Complex fast_pos_clogl(long double _Complex z)
{
   const long double rz = creall(z);
   const long double iz = cimagl(z);
   long double arg = atan2l(iz, rz);

   if (iz == 0.0L && arg < 0.0L)
      arg = -arg;

   return 0.5L*logl(rz*rz + iz*iz) + I*arg;
}
