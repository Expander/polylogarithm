!*********************************************************************
! This file is part of Polylogarithm.
!
! Polylogarithm is licenced under the MIT License.
!*********************************************************************

!*********************************************************************
!> @brief Clausen function \f$\operatorname{Cl}_4(\theta) = \operatorname{Im}(\operatorname{Li}_4(e^{i\theta}))\f$
!> @param x real angle
!> @return \f$\operatorname{Cl}_4(\theta)\f$
!> @author Alexander Voigt
!> Implemented as rational function approximation.
!*********************************************************************

double precision function dcl4(x)
  implicit none
  double precision :: x, y, z, z2, z4, p, q, p0, p1, h, sgn
  double precision, parameter :: PI = 3.14159265358979324D0
  double precision, parameter :: PI2 = 2*PI, PIH = PI/2, PI28 = PI*PI/8
  double precision, parameter :: zeta3 = 1.2020569031595943D0
  double precision, parameter :: cp(4) = (/ &
     -3.0555555555555556D-1,                &
      6.0521392328447206D-3,                &
     -1.9587493942041528D-5,                &
     -3.1137343767030358D-8                 /)
  double precision, parameter :: cq(4) = (/ &
      1.0000000000000000D+0,                &
     -2.2079728398400851D-2,                &
      1.0887447112236682D-4,                &
     -6.1847621370547954D-8                 /)
  double precision, parameter :: cr(6) = (/ &
      7.6223911686491336D-1,                &
     -2.4339587368267260D-1,                &
      2.8715364937979943D-2,                &
     -1.5368612510964667D-3,                &
      3.6261044225761673D-5,                &
     -2.8557977333851308D-7                 /)
  double precision, parameter :: cs(6) = (/ &
      1.0000000000000000D+0,                &
     -1.7465715261403233D-1,                &
      9.5439417991615653D-3,                &
     -1.7325070821666274D-4,                &
      5.9283675098376635D-7,                &
      9.4127575773361230D-10                /)

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
     dcl4 = 0
     return
  endif

  if (x .lt. PIH) then
    y = x*x
    z = y*y
    p = cp(1) + y * cp(2) + z * (cp(3) + y * cp(4))
    q = cq(1) + y * cq(2) + z * (cq(3) + y * cq(4))
    h = x*(zeta3 + y*(p/q + log(x)/6))
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

  dcl4 = sgn*h

end function dcl4
