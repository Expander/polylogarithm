! Compile as:
!
! gfortran example.F90 -I<install-directory>>/include <install-directory>/lib/libpolylogarithm_fortran.a

program example
  implicit none

#include "polylogarithm/Li2.fh"
#include "polylogarithm/Li3.fh"
#include "polylogarithm/Li4.fh"
#include "polylogarithm/Li5.fh"
#include "polylogarithm/Li6.fh"

  double precision :: x
  double complex :: z

  x = 1.1D0
  z = dcmplx(1.1D0, 1.1D0)

  print *, 'dli2(', x, ') = ', dli2(x)
  print *, 'dli3(', x, ') = ', dli3(x)
  print *, 'dli4(', x, ') = ', dli4(x)
  print *, 'cdli2(', z, ') = ', cdli2(z)
  print *, 'cdli3(', z, ') = ', cdli3(z)
  print *, 'cdli4(', z, ') = ', cdli4(z)
  print *, 'cdli5(', z, ') = ', cdli5(z)
  print *, 'cdli6(', z, ') = ', cdli6(z)

end program example
