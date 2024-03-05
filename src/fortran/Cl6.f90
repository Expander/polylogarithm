!*********************************************************************
! This file is part of Polylogarithm.
!
! Polylogarithm is licenced under the MIT License.
!*********************************************************************

!*********************************************************************
!> @brief Clausen function \f$\operatorname{Cl}_6(\theta) = \operatorname{Im}(\operatorname{Li}_6(e^{i\theta}))\f$
!> @param x real angle
!> @return \f$\operatorname{Cl}_6(\theta)\f$
!> @author Alexander Voigt
!> Implemented as rational function approximation.
!*********************************************************************

double precision function dcl6(x)
  implicit none
  double precision :: x, y, z, z2, z4, p, q, p0, p1, h, sgn
  double precision, parameter :: PI = 3.14159265358979324D0
  double precision, parameter :: PI2 = 2*PI, PIH = PI/2, PI28 = PI*PI/8
  double precision, parameter :: zeta3 = 1.2020569031595943D0
  double precision, parameter :: cp(4) = (/ &
      1.0369277551433699D+0,                &
     -2.0871954441071750D-1,                &
      2.0652251045312954D-2,                &
     -1.3834381382568400D-4                 /)
  double precision, parameter :: cq(4) = (/ &
      1.0000000000000000D+0,                &
     -8.0784096827362542D-3,                &
      5.8074568862993102D-6,                &
     -5.1960620033050114D-10                /)
  double precision, parameter :: cr(5) = (/ &
      7.9544504578027050D-1,                &
     -1.9255025309738589D-1,                &
      1.5805208288846591D-2,                &
     -5.4175380521534706D-4,                &
      6.7577493541009068D-6                 /)
  double precision, parameter :: cs(6) = (/ &
      1.0000000000000000D+0,                &
     -7.0798422394109274D-2,                &
      7.1744189715634762D-4,                &
      3.9098747334347093D-6,                &
      3.5669441618295266D-8,                &
      2.5315391843409925D-10                /)

  sgn = 1

  if (x .lt. 0) then
     x = -x
     sgn = -1
  endif

  if (x .ge. PI2) then
     x = mod(x, PI2)
  endif

  if (x .gt. PI) then
     p0 = 6.28125D0
     p1 = 0.0019353071795864769253D0
     x = (p0 - x) + p1
     sgn = -sgn
  endif

  if (x .eq. 0) then
     dcl6 = x
     return
  elseif (x .eq. PI) then
     dcl6 = 0
     return
  endif

  if (x .lt. PIH) then
    y = x*x
    z = y*y
    p = cp(1) + y * cp(2) + z * (cp(3) + y * cp(4))
    q = cq(1) + y * cq(2) + z * (cq(3) + y * cq(4))
    h = x*(p/q - 1.0D0/120*z*log(x))
  else
    y = PI - x
    z = y*y - PI28
    z2 = z*z
    z4 = z2*z2
    p = cr(1) + z * cr(2) + z2 * (cr(3) + z * cr(4)) + &
        z4 * cr(5)
    q = cs(1) + z * cs(2) + z2 * (cs(3) + z * cs(4)) + &
        z4 * (cs(5) + z * cs(6))
    h = y*p/q
  endif

  dcl6 = sgn*h

end function dcl6
