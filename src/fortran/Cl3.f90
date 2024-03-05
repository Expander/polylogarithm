!*********************************************************************
! This file is part of Polylogarithm.
!
! Polylogarithm is licenced under the MIT License.
!*********************************************************************

!*********************************************************************
!> @brief Clausen function \f$\operatorname{Cl}_3(\theta) = \operatorname{Re}(\operatorname{Li}_3(e^{i\theta}))\f$
!> @param x real angle
!> @return \f$\operatorname{Cl}_3(\theta)\f$
!> @author Alexander Voigt
!> Implemented as rational function approximation.
!*********************************************************************

double precision function dcl3(x)
  implicit none
  double precision :: x, y, z, z2, z4, p, q, p0, p1
  double precision, parameter :: PI = 3.14159265358979324D0
  double precision, parameter :: PI2 = 2*PI, PIH = PI/2, PI28 = PI*PI/8
  double precision, parameter :: zeta3 = 1.2020569031595943D0
  double precision, parameter :: cp(4) = (/ &
     -7.5000000000000001D-1,                &
      1.5707637881835541D-2,                &
     -3.5426736843494423D-5,                &
     -2.4408931585123682D-7                 /)
  double precision, parameter :: cq(4) = (/ &
      1.0000000000000000D+0,                &
     -2.5573146805410089D-2,                &
      1.5019774853075050D-4,                &
     -1.0648552418111624D-7                 /)
  double precision, parameter :: cr(6) = (/ &
     -4.9017024647634973D-1,                &
      4.1559155224660940D-1,                &
     -7.9425531417806701D-2,                &
      5.9420152260602943D-3,                &
     -1.8302227163540190D-4,                &
      1.8027408929418533D-6                 /)
  double precision, parameter :: cs(6) = (/ &
      1.0000000000000000D+0,                &
     -1.9495887541644712D-1,                &
      1.2059410236484074D-2,                &
     -2.5235889467301620D-4,                &
      1.0199322763377861D-6,                &
      1.9612106499469264D-9                 /)

  if (x .lt. 0) then
     x = -x
  endif

  if (x .ge. PI2) then
     x = mod(x, PI2)
  endif

  if (x .gt. PI) then
     p0 = 6.28125D0
     p1 = 0.0019353071795864769253D0
     x = (p0 - x) + p1
  endif

  if (x .eq. 0) then
     dcl3 = zeta3
  elseif (x .lt. PIH) then
     y = x*x
     z = y*y
     p = cp(1) + y * cp(2) + z * (cp(3) + y * cp(4))
     q = cq(1) + y * cq(2) + z * (cq(3) + y * cq(4))
     dcl3 = zeta3 + y*(p/q + log(x)/2)
  else
     y = PI - x
     z = y*y - PI28
     z2 = z*z
     z4 = z2*z2
     p = cr(1) + z * cr(2) + z2 * (cr(3) + z * cr(4)) + &
         z4 * (cr(5) + z * cr(6))
     q = cs(1) + z * cs(2) + z2 * (cs(3) + z * cs(4)) + &
         z4 * (cs(5) + z * cs(6))
     dcl3 = p/q
  endif

end function dcl3
