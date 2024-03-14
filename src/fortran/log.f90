!*********************************************************************
! This file is part of Polylogarithm.
!
! Polylogarithm is licenced under the MIT License.
!*********************************************************************


!*********************************************************************
!> @brief Implementation of log(1 + z) for complex z
!> @param z complex argument
!> @return log(1 + z)
!*********************************************************************
double complex function cdlog1p(z)
  implicit none
  double complex :: z, u
  double precision :: re, im

  u = 1 + z
  re = real(u)
  im = aimag(u)

  if (re .eq. 1 .and. im .eq. 0) then
     cdlog1p = z
  elseif (re .le. 0) then
     cdlog1p = log(u)
  else
     cdlog1p = log(u)*(z/(u - 1));
  endif

end function cdlog1p


!*********************************************************************
!> @brief Implementation of complex logarithm
!> @param z complex argument
!> @return log(z)
!> @note Points on the branch cut are treated differently from log(z):
!> Points with Im(z) == -0D0 are mapped to Im(z) == 0D0
!*********************************************************************
double complex function pos_cdlog(z)
  implicit none
  double complex :: z
  double precision :: re, im, arg

  re = real(z)
  im = aimag(z)

  if (im .eq. 0 .and. re .gt. 0) then
     pos_cdlog = dcmplx(log(re), 0.0D0)
  elseif (im .eq. 0) then
     pos_cdlog = dcmplx(log(-re), 3.14159265358979324D0)
  else
     pos_cdlog = log(z)
  endif

end function pos_cdlog
