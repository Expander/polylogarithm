Polylogarithm
=============

[![Build Status](https://travis-ci.org/Expander/polylogarithm.svg?branch=master)](https://travis-ci.org/Expander/polylogarithm)

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


Copying
-------

Polylogarithm is licenced under the GNU Lesser General Public License
(GNU LGPL) version 3.
