#include <math.h>

/**
 * @brief Clausen function \f$\operatorname{Cl}_2(\theta)\f$
 * @param x real angle
 * @return \f$\operatorname{Cl}_2(\theta)\f$
 * @author Alexander Voigt
 *
 * Implementation as series expansion in terms of Bernoulli numbers
 * [Abramowitz and Stegun, "Handbook of Mathematical Functions with
 * Formulas, Graphs, and Mathematical Tables", 27.8.2-3].
 */
double clausen_2_bernoulli(double x)
{
   const double PI = 3.14159265358979324;
   const double PI2 = 2*PI;
   const double A[] = {
      1.0000000000000000e+00, 1.3888888888888889e-02,
      6.9444444444444444e-05, 7.8735197782816830e-07,
      1.1482216343327454e-08, 1.8978869988970999e-10,
      3.3873013709535213e-12, 6.3726364431831804e-14,
      1.2462059912950672e-15, 2.5105444608999546e-17,
      5.1782588060906235e-19, 1.0887357368300849e-20,
      2.3257441143020872e-22, 5.0351952131473896e-24,
      1.1026499294381215e-25, 2.4386585509007345e-27,
      5.4401426788562523e-29, 1.2228340131217352e-30,
      2.7672634689679506e-32, 6.3000905918320139e-34,
      1.4420868388418475e-35, 3.3170939991595428e-37,
      7.6639135579206579e-39
   };
   const double B[] = {
       6.9314718055994531e-01, -4.1666666666666667e-02,
      -1.0416666666666667e-03, -4.9603174603174603e-05,
      -2.9279651675485009e-06, -1.9415383998717332e-07,
      -1.3870999114054670e-08, -1.0440290284867004e-09,
      -8.1670109639522231e-11, -6.5812165661369679e-12,
      -5.4297927275964755e-13, -4.5664875671936355e-14,
      -3.9019509040630692e-15, -3.3790622573736396e-16,
      -2.9599033551444005e-17, -2.6184896781186929e-18,
      -2.3365234885821292e-19, -2.1008128377954316e-20,
      -1.9016489757535847e-21, -1.7317557154340701e-22,
      -1.5855912475679039e-23, -1.4588733690004325e-24
   };

   double sgn = 1;

   if (x < 0) {
      x = -x;
      sgn = -1;
   }

   if (x >= PI2) {
      x = fmod(x, PI2);
   }

   if (x > PI) {
      const double p0 = 6.28125;
      const double p1 = 0.0019353071795864769253;
      x = (p0 - x) + p1;
      sgn = -sgn;
   }

   if (x == 0 || x == PI) {
      return 0;
   }

   double sum = 0;

   if (x < PI/2) {
      const double x2 = x*x;
      sum = A[sizeof(A)/sizeof(A[0]) - 1];
      for (int i = sizeof(A)/sizeof(A[0]) - 2; i >= 0; --i) {
         sum = x2*sum + A[i];
      }
      sum -= log(x);
      sum *= x;
   } else {
      const double d = PI - x;
      const double d2 = d*d;
      sum = B[sizeof(B)/sizeof(B[0]) - 1];
      for (int i = sizeof(B)/sizeof(B[0]) - 2; i >= 0; --i) {
         sum = d2*sum + B[i];
      }
      sum *= d;
   }

   return sgn*sum;
}
