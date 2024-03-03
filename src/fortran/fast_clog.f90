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
  fast_cdlog = dcmplx(log(hypot(re, im)), datan2(im, re))

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

  if (im .eq. 0 .and. re .gt. 0) then
     fast_pos_cdlog = dcmplx(log(re), 0.0D0)
  elseif (im .eq. 0) then
     fast_pos_cdlog = dcmplx(log(-re), 3.14159265358979324D0)
  else
     fast_pos_cdlog = dcmplx(log(hypot(re, im)), datan2(im, re))
  endif

end function fast_pos_cdlog
