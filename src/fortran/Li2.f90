!*********************************************************************
! This file is part of Polylogarithm.
!
! Polylogarithm is licenced under the MIT License.
!*********************************************************************


!*********************************************************************
!> @brief Real dilogarithm \f$\operatorname{Li}_2(x)\f$
!> @param x real argument
!> @return \f$\operatorname{Li}_2(x)\f$
!> @author Alexander Voigt
!>
!> Implemented as a rational function approximation with a maximum
!> error of 5e-17
!> [[arXiv:2201.01678](https://arxiv.org/abs/2201.01678)].
!*********************************************************************

double precision function dli2(x)
  implicit none
  double precision :: x, y, r, s, y2, y4, p, q, l
  double precision, parameter :: PI = 3.14159265358979324D0
  double precision, parameter :: cp(6) = (/ &
      0.9999999999999999502D+0,             &
     -2.6883926818565423430D+0,             &
      2.6477222699473109692D+0,             &
     -1.1538559607887416355D+0,             &
      2.0886077795020607837D-1,             &
     -1.0859777134152463084D-2              /)
  double precision, parameter :: cq(7) = (/ &
      1.0000000000000000000D+0,             &
     -2.9383926818565635485D+0,             &
      3.2712093293018635389D+0,             &
     -1.7076702173954289421D+0,             &
      4.1596017228400603836D-1,             &
     -3.9801343754084482956D-2,             &
      8.2743668974466659035D-4              /)

  ! transform to [0, 1/2]
   if (x .lt. -1) then
      l = log(1 - x)
      y = 1/(1 - x)
      r = -PI**2/6 + l*(0.5D0*l - log(-x))
      s = 1
   elseif (x .eq. -1) then
      dli2 = -PI**2/12
      return
   elseif (x .lt. 0) then
      y = x/(x - 1)
      r = -0.5D0*log(1 - x)**2
      s = -1
   elseif (x .eq. 0) then
      dli2 = x
      return
   elseif (x .lt. 0.5D0) then
      y = x
      r = 0
      s = 1
   elseif (x .lt. 1) then
      y = 1 - x
      r = PI**2/6 - log(x)*log(y)
      s = -1
   elseif (x .eq. 1) then
      dli2 = PI**2/6
      return
   elseif (x .lt. 2) then
      l = log(x)
      y = 1 - 1/x
      r = PI**2/6 - l*(log(y) + 0.5D0*l)
      s = 1
   else
      y = 1/x
      r = PI**2/3 - 0.5D0*log(x)**2
      s = -1
   endif

  y2 = y*y
  y4 = y2*y2
  p = cp(1) + y * cp(2) + y2 * (cp(3) + y * cp(4)) +      &
      y4 * (cp(5) + y * cp(6))
  q = cq(1) + y * cq(2) + y2 * (cq(3) + y * cq(4)) +      &
      y4 * (cq(5) + y * cq(6) + y2 * cq(7))

  dli2 = r + s*y*p/q

end function dli2


!*********************************************************************
!> @brief Complex dilogarithm \f$\operatorname{Li}_2(z)\f$
!> @param z complex argument
!> @return \f$\operatorname{Li}_2(z)\f$
!> @note Implementation translated from SPheno by Alexander Voigt
!*********************************************************************

double complex function cdli2(z)
  implicit none
  double complex :: z, rest, u, u2, u4, sum, cdlog1p
  double precision :: rz, iz, nz, sgn, dli2
  double precision, parameter :: PI = 3.14159265358979324D0
  double precision, parameter :: bf(10) = (/ &
    - 1.0D0/4.0D0,                           &
    + 1.0D0/36.0D0,                          &
    - 1.0D0/3600.0D0,                        &
    + 1.0D0/211680.0D0,                      &
    - 1.0D0/10886400.0D0,                    &
    + 1.0D0/526901760.0D0,                   &
    - 4.0647616451442255D-11,                &
    + 8.9216910204564526D-13,                &
    - 1.9939295860721076D-14,                &
    + 4.5189800296199182D-16                 /)

  rz = real(z)
  iz = aimag(z)

  ! special cases
  if (iz .eq. 0) then
     if (rz .le. 1) cdli2 = dcmplx(dli2(rz), iz)
     if (rz .gt. 1) cdli2 = dcmplx(dli2(rz), -PI*log(rz))
     return
  endif

  nz = rz**2 + iz**2

  if (nz .lt. EPSILON(1D0)) then
     cdli2 = z*(1 + 0.25D0*z)
     return
  endif

  ! transformation to |z| < 1, Re(z) <= 0.5
  if (rz .le. 0.5D0) then
     if (nz .gt. 1) then
        u = -cdlog1p(-1/z)
        rest = -0.5D0*log(-z)**2 - PI**2/6
        sgn = -1
     else ! nz <= 1
        u = -cdlog1p(-z)
        rest = 0
        sgn = 1
     endif
  else ! rz > 0.5D0
     if (nz .le. 2*rz) then
        u = -log(z)
        rest = u*cdlog1p(-z) + PI**2/6
        sgn = -1
     else ! nz > 2*rz
        u = -cdlog1p(-1/z)
        rest = -0.5D0*log(-z)**2 - PI**2/6
        sgn = -1
     endif
  endif

  u2 = u**2
  u4 = u2**2
  sum =                                                    &
     u +                                                   &
     u2 * (bf(1) +                                         &
     u  * (bf(2) +                                         &
     u2 * (                                                &
         bf(3) +                                           &
         u2*bf(4) +                                        &
         u4*(bf(5) + u2*bf(6)) +                           &
         u4*u4*(bf(7) + u2*bf(8) + u4*(bf(9) + u2*bf(10))) &
     )))

  cdli2 = sgn*sum + rest

end function cdli2
