!*********************************************************************
! This file is part of Polylogarithm.
!
! Polylogarithm is licenced under the MIT License.
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
  fast_cdlog = dcmplx(0.5D0*log(re**2 + im**2), datan2(im, re))

end function fast_cdlog


!*********************************************************************
!> @brief Fast implementation of complex logarithm
!> @param z complex argument
!> @return log(z)
!> @note Points on the branch cut are treated differently from log(z):
!> Points with Im(z) == -0D0 are mapped to Im(z) == 0D0
!*********************************************************************
double complex function fast_pos_cdlog(z)
  implicit none
  double complex :: z
  double precision :: re, im, arg

  re = real(z)
  im = aimag(z)
  arg = datan2(im, re)

  if (im .eq. 0 .and. arg .lt. 0) arg = -arg

  fast_pos_cdlog = dcmplx(0.5D0*log(re**2 + im**2), arg)

end function fast_pos_cdlog
