/*
  Compile as follows:

  gcc example.c <install-directory>/lib/libpolylogarithm_c.a -lm
 */

#include <complex.h>
#include <stdio.h>

double li2(double x);
double _Complex cli2(double _Complex z);

int main() {
   const double x = 1.1;
   printf("li2(%g) = %g\n", x, li2(x));

   const double _Complex z = 1.1 + 1.1*I;
   printf("li2(%g + %gi) = %g + %gi\n", creal(z), cimag(z), creal(cli2(z)), cimag(cli2(z)));
}
