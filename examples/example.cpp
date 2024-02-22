// Compile as follows:
//
// g++ examples/example.cpp -I<install-directory>/include <install-directory>/lib/libpolylogarithm_cpp.a


#include <complex>
#include <iostream>
#include "polylogarithm/Cl.hpp"
#include "polylogarithm/Cl1.hpp"
#include "polylogarithm/Cl2.hpp"
#include "polylogarithm/Cl3.hpp"
#include "polylogarithm/Cl4.hpp"
#include "polylogarithm/Cl5.hpp"
#include "polylogarithm/Cl6.hpp"
#include "polylogarithm/Li.hpp"
#include "polylogarithm/Li2.hpp"
#include "polylogarithm/Li3.hpp"
#include "polylogarithm/Li4.hpp"
#include "polylogarithm/Li5.hpp"
#include "polylogarithm/Li6.hpp"
#include "polylogarithm/Sl.hpp"

int main() {
   const double x = 1.1;
   std::cout << "Cl1(" << x << ") = " << polylogarithm::Cl1(x) << '\n';
   std::cout << "Cl2(" << x << ") = " << polylogarithm::Cl2(x) << '\n';
   std::cout << "Cl3(" << x << ") = " << polylogarithm::Cl3(x) << '\n';
   std::cout << "Cl4(" << x << ") = " << polylogarithm::Cl4(x) << '\n';
   std::cout << "Cl5(" << x << ") = " << polylogarithm::Cl5(x) << '\n';
   std::cout << "Cl6(" << x << ") = " << polylogarithm::Cl6(x) << '\n';
   std::cout << "Cl(10," << x << ") = " << polylogarithm::Cl(10,x) << '\n';
   std::cout << "Li2(" << x << ") = " << polylogarithm::Li2(x) << '\n';
   std::cout << "Li3(" << x << ") = " << polylogarithm::Li3(x) << '\n';
   std::cout << "Li4(" << x << ") = " << polylogarithm::Li4(x) << '\n';
   std::cout << "Li(10," << x << ") = " << polylogarithm::Li(10,x) << '\n';
   std::cout << "Sl(10," << x << ") = " << polylogarithm::Sl(10,x) << '\n';

   const std::complex<double> z(1.1, 1.1);
   std::cout << "Li2(" << z << ") = " << polylogarithm::Li2(z) << '\n';
   std::cout << "Li3(" << z << ") = " << polylogarithm::Li3(z) << '\n';
   std::cout << "Li4(" << z << ") = " << polylogarithm::Li4(z) << '\n';
   std::cout << "Li5(" << z << ") = " << polylogarithm::Li5(z) << '\n';
   std::cout << "Li6(" << z << ") = " << polylogarithm::Li6(z) << '\n';
   std::cout << "Li(10," << z << ") = " << polylogarithm::Li(10,z) << '\n';
}
