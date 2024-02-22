/*
  Compile as follows:

  gcc example.c -I<install-directory>/include <install-directory>/lib/libpolylogarithm_c.a -lm
 */

#include <complex.h>
#include <stdio.h>
#include "polylogarithm/Cl2.h"
#include "polylogarithm/Cl3.h"
#include "polylogarithm/Cl4.h"
#include "polylogarithm/Cl5.h"
#include "polylogarithm/Cl6.h"
#include "polylogarithm/Li2.h"
#include "polylogarithm/Li3.h"
#include "polylogarithm/Li4.h"
#include "polylogarithm/Li5.h"
#include "polylogarithm/Li6.h"

int main() {
   const double x = 1.1;
   printf("cl2(%g) = %g\n", x, cl2(x));
   printf("cl3(%g) = %g\n", x, cl3(x));
   printf("cl4(%g) = %g\n", x, cl4(x));
   printf("cl5(%g) = %g\n", x, cl5(x));
   printf("cl6(%g) = %g\n", x, cl6(x));
   printf("li2(%g) = %g\n", x, li2(x));
   printf("li3(%g) = %g\n", x, li3(x));
   printf("li4(%g) = %g\n", x, li4(x));

   const double _Complex z = 1.1 + 1.1*I;
   printf("cli2(%g + %gi) = %g + %gi\n", creal(z), cimag(z), creal(cli2(z)), cimag(cli2(z)));
   printf("cli3(%g + %gi) = %g + %gi\n", creal(z), cimag(z), creal(cli3(z)), cimag(cli3(z)));
   printf("cli4(%g + %gi) = %g + %gi\n", creal(z), cimag(z), creal(cli4(z)), cimag(cli4(z)));
   printf("cli5(%g + %gi) = %g + %gi\n", creal(z), cimag(z), creal(cli5(z)), cimag(cli5(z)));
   printf("cli6(%g + %gi) = %g + %gi\n", creal(z), cimag(z), creal(cli6(z)), cimag(cli6(z)));
}
