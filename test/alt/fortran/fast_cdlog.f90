!*********************************************************************
! This file is part of Polylogarithm.
!
! Polylogarithm is licenced under the GNU Lesser General Public
! License (GNU LGPL) version 3.
!*********************************************************************


!*********************************************************************
!> @brief Fast implementation of complex logarithm
!> @param z complex argument
!> @return log(z)
!*********************************************************************
double complex function fast_cdlog(z)
  implicit none
  double complex :: z
  double precision :: re, im

  re = real(z)
  im = aimag(z)
  fast_cdlog = dcmplx(0.5D0*log(re*re + im*im), datan2(im, re))

end function fast_cdlog
