/* ====================================================================
 * This file is part of Polylogarithm.
 *
 * Polylogarithm is licenced under the MIT License.
 * ==================================================================== */

#include <complex.h>
#include <math.h>


double _Complex clog1p(double _Complex z)
{
   const double _Complex u = 1.0 + z;
   const double rz = creal(u);
   const double iz = cimag(u);

   if (rz == 1.0 && iz == 0.0) {
      return z;
   } else if (rz <= 0.0) {
      return clog(u);
   }

   return clog(u)*(z/(u - 1.0));
}


long double _Complex clog1pl(long double _Complex z)
{
   const long double _Complex u = 1.0L + z;
   const long double rz = creall(u);
   const long double iz = cimagl(u);

   if (rz == 1.0L && iz == 0.0L) {
      return z;
   } else if (rz <= 0.0L) {
      return clogl(u);
   }

   return clogl(u)*(z/(u - 1.0L));
}


double _Complex pos_clog(double _Complex z)
{
   const double rz = creal(z);
   const double iz = cimag(z);

   if (iz == 0.0 && rz > 0.0) {
      return log(rz);
   } else if (iz == 0.0) {
      return log(-rz) + I*3.1415926535897932;
   }

   return clog(z);
}


long double _Complex pos_clogl(long double _Complex z)
{
   const long double rz = creall(z);
   const long double iz = cimagl(z);

   if (iz == 0.0L && rz > 0.0L) {
      return logl(rz);
   } else if (iz == 0.0L) {
      return logl(-rz) + I*3.14159265358979323846264338327950288L;
   }

   return clogl(z);
}
