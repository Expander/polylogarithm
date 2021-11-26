Polylogarithm
=============

[![Build Status](https://github.com/Expander/polylogarithm/workflows/test/badge.svg)](https://github.com/Expander/polylogarithm/actions)

The Polylogarithm package provides C, C++ and Fortran implementations
of various polylogarithms, including the real and complex dilogarithm,
trilogarithm, and standard Clausen functions.  The implementations
have been fully tested against the literature and many other
implementations and are highly optimized for fast numerical
evaluation.

The package has no external dependencies, except for the C/C++/Fortran
standard libraries.  The implementations of the individual polylogarithm
functions are distributed among different source code files, so
individual source code files can be easily extracted and incorporated
into existing projects.


Notes
-----

The implementation of the complex dilogarithm is inspired by the
implementation in [SPheno](https://spheno.hepforge.org/).


Citation
--------

~~~
@software{polylogarithm,
    author       = {{Alexander Voigt}},
    title        = {{Polylogarithm}},
    year         = {2021},
    version      = {6.6.0},
    url          = {https://github.com/Expander/polylogarithm},
    note         = {[License: MIT]}
}
~~~


Copying
-------

Polylogarithm is licenced under the MIT License.
