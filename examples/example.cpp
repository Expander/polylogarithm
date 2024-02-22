// Compile as follows:
//
// g++ examples/example.cpp -I<install-directory>/include <install-directory>/lib/libpolylogarithm_cpp.a


#include <complex>
#include <iostream>
#include "polylogarithm/Li2.hpp"

int main() {
   const double x = 1.1;
   std::cout << "Li2(" << x << ") = " << polylogarithm::Li2(x) << '\n';

   const std::complex<double> z(1.1, 1.1);
   std::cout << "Li2(" << z << ") = " << polylogarithm::Li2(z) << '\n';
}
