! Compile as:
!
! gfortran example.f90 <install-directory>/lib/libpolylogarithm_fortran.a

program example
  implicit none

  double precision :: x, dli2
  double complex :: z, cdli2

  x = 1.1D0
  z = dcmplx(1.1D0, 1.1D0)

  print *, 'dli2(', x, ') = ', dli2(x)
  print *, 'cdli2(', z, ') = ', cdli2(z)

end program example
