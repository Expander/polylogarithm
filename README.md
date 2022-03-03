Polylogarithm
=============

[![Build Status](https://github.com/Expander/polylogarithm/workflows/test/badge.svg)](https://github.com/Expander/polylogarithm/actions)

The Polylogarithm package provides C, C++ and Fortran implementations
of various polylogarithms, including the real and complex dilogarithm,
trilogarithm, and (Standard and Glaisher) Clausen functions.  The
implementations have been fully tested against the literature and many
other implementations and are highly optimized for fast numerical
evaluation.

The package has no external dependencies, except for the C/C++/Fortran
standard libraries.  The implementations of the individual polylogarithm
functions are distributed among different source code files, so
individual source code files can be easily extracted and incorporated
into existing projects.


Example in C++
--------------

```.cpp
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
```

Output:

```
Li_2(1) = 1.64493
Li_3(1) = 1.20206
Li_4(1) = 1.08232
Li_2((1,1)) = (0.61685,1.46036)
Li_3((1,1)) = (0.871159,1.26708)
Li_4((1,1)) = (0.959319,1.13804)
Li_5((1,1)) = (0.987467,1.06844)
Li_6((1,1)) = (0.99615,1.03355)
Li_10((1,1)) = (0.999962,1.00199)
```


Notes
-----

The implementation of the real dilogarithm is an adaption of
[[arXiv:2201.01678](https://arxiv.org/abs/2201.01678)].

The implementation of the complex dilogarithm is inspired by the
implementation in [SPheno](https://spheno.hepforge.org/).

The implementation of the general n-th order complex polylogarithm an
adaption of [[arXiv:2010.09860](https://arxiv.org/abs/2010.09860)].

Citation
--------

```.bibtex
@software{polylogarithm,
    author       = {{Alexander Voigt}},
    title        = {{Polylogarithm}},
    year         = {2022},
    version      = {6.11.0},
    url          = {https://github.com/Expander/polylogarithm},
    note         = {[License: MIT]}
}
```


Copying
-------

Polylogarithm is licenced under the MIT License.
