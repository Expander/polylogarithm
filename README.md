Polylogarithm
=============

[![Build Status](https://travis-ci.org/Expander/polylogarithm.svg?branch=master)](https://travis-ci.org/Expander/polylogarithm)

The Polylogarithm package provides C++ implementations of various
polylogarithms, including the real and complex dilogarithm,
trilogarithm, and standard Clausen functions.

The package has no external dependencies, except for the C++ standard
library.  The implementations of the individual polylogarithm
functions are distributed among different source code files and are
completely separate from each other, so individual source code files
can be easily extracted and incorporated into existing projects.


Notes
-----

The implementation of the real dilogarithm (with long double
precision) is inspired by the implementation from the ROOT package
(root.cern.ch), licensed under the GNU LGPL, which is based on a
FORTRAN implementation by K.S. Kölbig (CERNLIB DILOG function C332).

The implementation of the complex dilogarithm is inspired by the
implementation in SPheno (spheno.hepforge.org).


Copying
-------

Polylogarithm is licenced under the GNU Lesser General Public License
(GNU LGPL) version 3.
