#include <complex.h>
#include <math.h>


static long double complex Li2series(long double complex x1)
{
   /* these are the even-n Bernoulli numbers, already divided by (n + 1)!
      as in Table[BernoulliB[n]/(n + 1)!, {n, 2, 50, 2}] */
   long double b[25] = {
       0.02777777777777777777777777777777777777777777778e0l,
      -0.000277777777777777777777777777777777777777777778e0l,
       4.72411186696900982615268329554043839758125472e-6l,
      -9.18577307466196355085243974132863021751910641e-8l,
       1.89788699889709990720091730192740293750394761e-9l,
      -4.06476164514422552680590938629196667454705711e-11l,
       8.92169102045645255521798731675274885151428361e-13l,
      -1.993929586072107568723644347793789705630694749e-14l,
       4.51898002961991819165047655285559322839681901e-16l,
      -1.035651761218124701448341154221865666596091238e-17l,
       2.39521862102618674574028374300098038167894899e-19l,
      -5.58178587432500933628307450562541990556705462e-21l,
       1.309150755418321285812307399186592301749849833e-22l,
      -3.087419802426740293242279764866462431595565203e-24l,
       7.31597565270220342035790560925214859103339899e-26l,
      -1.740845657234000740989055147759702545340841422e-27l,
       4.15763564461389971961789962077522667348825413e-29l,
      -9.96214848828462210319400670245583884985485196e-31l,
       2.394034424896165300521167987893749562934279156e-32l,
      -5.76834735536739008429179316187765424407233225e-34l,
       1.393179479647007977827886603911548331732410612e-35l,
      -3.372121965485089470468473635254930958979742891e-37l,
       8.17820877756210262176477721487283426787618937e-39l,
      -1.987010831152385925564820669234786567541858996e-40l,
       4.83577851804055089628705937311537820769430091e-42l
   };

   long double complex xm = -clogl(x1);
   long double complex x2 = xm*xm;
   long double complex ls = xm - x2/4;

   for (int j = 0; j < 25; j++) {
      long double complex ls0 = ls;
      xm = xm*x2;
      ls = ls + xm*b[j];
      long double complex ls1 = ls;
      if (ls0 == ls1)
         break;
   }

   return ls;
}


/*
  Implementation of dilogarithm from FeynHiggs 2.16.0.

  file: src/LT/spence.F

  Translated to C by Alexander Voigt
 */
static long double complex spence(int i_in, long double complex x_in, int s)
{
   long double complex x[2], sp;
   long double I1 = 1;
   long double I2 = I1/2;
   long double pi = 3.1415926535897932384626433832795029l;
   long double zeta2 = pi*pi/6;
   long double zeroeps = 1e-20l;
   long double complex cI = 0.0l + 1.0lI;
   long double complex cIeps = cI*1e-50l;

   x[i_in] = x_in;
   x[1 - i_in] = 1 - x_in;

   if (creall(x[0]) < I2) {
      if (cabsl(x[0]) < 1) {
         sp = Li2series(x[1] - s * cIeps);
      } else {
         long double complex l = clogl(-x[0] - s * cIeps);
         sp = -zeta2 - l * l / 2 - Li2series(-x[1] / x[0] + s * cIeps);
      }
   } else {
      sp = zeta2;
      long double ax1 = cabsl(x[1]);

      if (ax1 > zeroeps) {
         sp = zeta2 - clogl(x[0] + s * cIeps) * clogl(x[1] - s * cIeps);
         if (ax1 < 1) {
            sp = sp - Li2series(x[0] + s * cIeps);
         } else {
            long double complex l = clogl(-x[1] - s * cIeps);
            sp = sp + zeta2 + l * l / 2 + Li2series(-x[0] / x[1] - s * cIeps);
         }
      }
   }

   return sp;
}


void looptools_dilog(long double re, long double im, long double* res_re, long double* res_im)
{
   int zPrec = 0;
   long double complex z = re + I*im;
   long double complex result = spence(0, z, zPrec);
   *res_re = creall(result);
   *res_im = cimagl(result);
}
