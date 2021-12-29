!*********************************************************************
! This file is part of Polylogarithm.
!
! Polylogarithm is licenced under the MIT License.
!*********************************************************************

!*********************************************************************
!> @brief Clausen function \f$\operatorname{Cl}_2(\theta) = \operatorname{Im}(\operatorname{Li}_2(e^{i\theta}))\f$
!> @param x real angle
!> @return \f$\operatorname{Cl}_2(\theta)\f$
!> @author Alexander Voigt
!> Implemented as economized Padé approximation.
!*********************************************************************

double precision function dcl2(x)
  implicit none
  double precision :: x, y, z, z2, z4, p, q, p0, p1, h, sgn
  double precision, parameter :: PI = 3.14159265358979324D0
  double precision, parameter :: PI2 = 2*PI, PIH = PI/2, PI28 = PI*PI/8
  double precision, parameter :: cp(4) = (/ &
      1.3888888888888889D-2,                &
     -4.3286930203743071D-4,                &
      3.2779814789973427D-6,                &
     -3.6001540369575084D-9                 /)
  double precision, parameter :: cq(4) = (/ &
      1.0000000000000000D+0,                &
     -3.6166589746694121D-2,                &
      3.6015827281202639D-4,                &
     -8.3646182842184428D-7                 /)
  double precision, parameter :: cr(6) = (/ &
      6.4005702446195512D-1,                &
     -2.0641655351338783D-1,                &
      2.4175305223497718D-2,                &
     -1.2355955287855728D-3,                &
      2.5649833551291124D-5,                &
     -1.4783829128773320D-7                 /)
  double precision, parameter :: cs(6) = (/ &
      1.0000000000000000D+0,                &
     -2.5299102015666356D-1,                &
      2.2148751048467057D-2,                &
     -7.8183920462457496D-4,                &
      9.5432542196310670D-6,                &
     -1.8184302880448247D-8                 /)

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

  if (x .eq. 0 .or. x .eq. PI) then
     dcl2 = 0
     return
  endif

  if (x .lt. PIH) then
    y = x*x
    z = y*y
    p = cp(1) + y * cp(2) + z * (cp(3) + y * cp(4))
    q = cq(1) + y * cq(2) + z * (cq(3) + y * cq(4))
    h = x*(1 - log(x) + y*p/q)
  else
    y = PI - x
    z = y*y - PI28
    z2 = z*z
    z4 = z2*z2
    p = cr(1) + z * cr(2) + z2 * (cr(3) + z * cr(4)) + &
        z4 * (cr(5) + z * cr(6))
    q = cs(1) + z * cs(2) + z2 * (cs(3) + z * cs(4)) + &
        z4 * (cs(5) + z * cs(6))
    h = y*p/q
  endif

  dcl2 = sgn*h

end function dcl2
