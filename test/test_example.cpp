#include "Li.hpp"
#include "Li2.hpp"
#include "Li3.hpp"
#include "Li4.hpp"
#include "Li5.hpp"
#include "Li6.hpp"
#include <iostream>

int main() {
   using namespace polylogarithm;

   const double x = 1.0;
   const std::complex<double> z(1.0, 1.0);

   // real polylogarithms for real arguments
   std::cout
      << "Li_2(" << x << ") = " << Li2(x) << '\n'
      << "Li_3(" << x << ") = " << Li3(x) << '\n'
      << "Li_4(" << x << ") = " << Li4(x) << '\n';

   // complex polylogarithms for complex arguments
   std::cout
      << "Li_2(" << z << ") = " << Li2(z) << '\n'
      << "Li_3(" << z << ") = " << Li3(z) << '\n'
      << "Li_4(" << z << ") = " << Li4(z) << '\n'
      << "Li_5(" << z << ") = " << Li5(z) << '\n'
      << "Li_6(" << z << ") = " << Li6(z) << '\n'
      << "Li_10(" << z << ") = " << Li(10,z) << '\n';
}
