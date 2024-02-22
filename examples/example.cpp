// Compile as follows:
//
// g++ examples/example.cpp -I<install-directory>/include <install-directory>/lib/libpolylogarithm_cpp.a


#include <complex>
#include <iostream>
#include "polylogarithm/Li2.hpp"
#include "polylogarithm/Li3.hpp"
#include "polylogarithm/Li4.hpp"
#include "polylogarithm/Li5.hpp"
#include "polylogarithm/Li6.hpp"

int main() {
   const double x = 1.1;
   std::cout << "Li2(" << x << ") = " << polylogarithm::Li2(x) << '\n';
   std::cout << "Li3(" << x << ") = " << polylogarithm::Li3(x) << '\n';
   std::cout << "Li4(" << x << ") = " << polylogarithm::Li4(x) << '\n';

   const std::complex<double> z(1.1, 1.1);
   std::cout << "Li2(" << z << ") = " << polylogarithm::Li2(z) << '\n';
   std::cout << "Li3(" << z << ") = " << polylogarithm::Li3(z) << '\n';
   std::cout << "Li4(" << z << ") = " << polylogarithm::Li4(z) << '\n';
   std::cout << "Li5(" << z << ") = " << polylogarithm::Li5(z) << '\n';
   std::cout << "Li6(" << z << ") = " << polylogarithm::Li6(z) << '\n';
}
