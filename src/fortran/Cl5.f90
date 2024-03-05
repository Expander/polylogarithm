!*********************************************************************
! This file is part of Polylogarithm.
!
! Polylogarithm is licenced under the MIT License.
!*********************************************************************

!*********************************************************************
!> @brief Clausen function \f$\operatorname{Cl}_5(\theta) = \operatorname{Re}(\operatorname{Li}_5(e^{i\theta}))\f$
!> @param x real angle
!> @return \f$\operatorname{Cl}_5(\theta)\f$
!> @author Alexander Voigt
!> Implemented as a rational function approximation.
!*********************************************************************

double precision function dcl5(x)
  implicit none
  double precision :: x, y, z, z2, z4, p, q, p0, p1
  double precision, parameter :: PI = 3.14159265358979324D0
  double precision, parameter :: PI2 = 2*PI, PIH = PI/2, PI28 = PI*PI/8
  double precision, parameter :: zeta5 = 1.0369277551433699D0
  double precision, parameter :: cp(4) = (/ &
      1.0369277551433699D+0,                &
     -6.1354800479984468D-1,                &
      9.4076401395712763D-2,                &
     -9.4056155866704436D-4                 /)
  double precision, parameter :: cq(5) = (/ &
      1.0000000000000000D+0,                &
     -1.2073698633244778D-2,                &
      1.3703409625482991D-5,                &
     -1.9701280330628469D-9,                &
      2.1944550184416500D-11                /)
  double precision, parameter :: cr(6) = (/ &
     -4.5930112735784898D-1,                &
      4.3720705508867954D-1,                &
     -7.5895226486465095D-2,                &
      5.2244176912488065D-3,                &
     -1.5677716622013956D-4,                &
      1.6641624171748576D-6                 /)
  double precision, parameter :: cs(6) = (/ &
      1.0000000000000000D+0,                &
     -1.2211486825401188D-1,                &
      3.8940070749313620D-3,                &
     -2.2674805547074318D-5,                &
     -7.4383354448335299D-8,                &
     -3.4131758392216437D-10                /)

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
     dcl5 = zeta5
  elseif (x .lt. PIH) then
     y = x*x
     z = y*y
     p = cp(1) + y * cp(2) + z * (cp(3) + y * cp(4))
     q = cq(1) + y * cq(2) + z * (cq(3) + y * cq(4) + z * cq(5))
     dcl5 = p/q - 1.0D0/24*z*log(x);
  else
     y = PI - x
     z = y*y - PI28
     z2 = z*z
     z4 = z2*z2
     p = cr(1) + z * cr(2) + z2 * (cr(3) + z * cr(4)) + &
         z4 * (cr(5) + z * cr(6))
     q = cs(1) + z * cs(2) + z2 * (cs(3) + z * cs(4)) + &
         z4 * (cs(5) + z * cs(6))
     dcl5 = p/q
  endif

end function dcl5
