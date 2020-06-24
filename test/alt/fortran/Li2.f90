!******************************************************************************
!> @brief Evaluation of polynomial P(x) with len coefficients c
!> @param x real argument of P
!> @param c coefficients of P(x)
!> @param len number of coefficients
!> @return \f$\mathrm{Li}_2(x)\f$
!******************************************************************************

double precision function dhorner(x, c, len)
  implicit none
  integer :: len, i
  double precision :: x
  double precision, dimension(len) :: c

  dhorner = 0

  do i = len, 1, -1
     dhorner = dhorner*x + c(i)
  end do

end function dhorner


!******************************************************************************
!> @brief Fast implementation of complex logarithm
!> @param z complex argument
!> @return log(z)
!******************************************************************************
double complex function fast_cdlog(z)
  implicit none
  double complex :: z
  double precision :: re, im

  re = real(z)
  im = aimag(z)
  fast_cdlog = dcmplx(0.5D0*log(re*re + im*im), datan2(im, re))

end function fast_cdlog


!******************************************************************************
!> @brief Real dilogarithm \f$\mathrm{Li}_2(x)\f$
!> @param x real argument
!> @return \f$\mathrm{Li}_2(x)\f$
!> @author Alexander Voigt
!>
!> Implemented as an economized Pade approximation with a
!> maximum error of 4.16e-18.
!******************************************************************************

double precision function dli2(x)
  implicit none
  double precision, parameter :: PI = 3.14159265358979324D0
  double precision :: x, y, r, s, z, p, q, l, dhorner
  double precision, dimension(8) :: cp, cq

  cp(1) =  1.0706105563309304277D+0
  cp(2) = -4.5353562730201404017D+0
  cp(3) =  7.4819657596286408905D+0
  cp(4) = -6.0516124315132409155D+0
  cp(5) =  2.4733515209909815443D+0
  cp(6) = -4.6937565143754629578D-1
  cp(7) =  3.1608910440687221695D-2
  cp(8) = -2.4630612614645039828D-4

  cq(1) =  1.0000000000000000000D+0
  cq(2) = -4.5355682121856044935D+0
  cq(3) =  8.1790029773247428573D+0
  cq(4) = -7.4634190853767468810D+0
  cq(5) =  3.6245392503925290187D+0
  cq(6) = -8.9936784740041174897D-1
  cq(7) =  9.8554565816757007266D-2
  cq(8) = -3.2116618742475189569D-3

  ! transform to [0, 1/2)
   if (x .lt. -1) then
      l = log(1 - x)
      y = 1/(1 - x)
      r = -PI*PI/6 + l*(0.5D0*l - log(-x))
      s = 1
   elseif (x .eq. -1) then
      dli2 = -PI*PI/12
      return
   elseif (x .lt. 0) then
      y = x/(x - 1)
      r = -0.5D0*log(1 - x)**2
      s = -1
   elseif (x .eq. 0) then
      dli2 = 0
      return
   elseif (x .lt. 0.5D0) then
      y = x
      r = 0
      s = 1
   elseif (x .lt. 1) then
      y = 1 - x
      r = PI*PI/6 - log(x)*log(1 - x)
      s = -1
   elseif (x .eq. 1) then
      dli2 = PI*PI/6
      return
   elseif (x .lt. 2) then
      l = log(x)
      y = 1 - 1/x
      r = PI*PI/6 - l*(log(1 - 1/x) + 0.5D0*l)
      s = 1
   else
      y = 1/x
      r = PI*PI/3 - 0.5D0*log(x)**2
      s = -1
   endif

  z = y - 0.25D0

  p = dhorner(z, cp, 8)
  q = dhorner(z, cq, 8)

  dli2 = r + s*y*p/q

end function dli2


!******************************************************************************
!> @brief Complex dilogarithm \f$\mathrm{Li}_2(z)\f$
!> @param z complex argument
!> @note Implementation translated from SPheno by Alexander Voigt
!> @return \f$\mathrm{Li}_2(z)\f$
!******************************************************************************
double complex function cdli2(z)
  implicit none
  double precision, parameter :: PI = 3.14159265358979324D0
  double precision :: rz, iz, nz, sgn, dli2
  double complex :: z, cy, cz, cz2, sum, fast_cdlog
  double precision, dimension(10) :: bf

  bf( 1) = - 1.0D0/4.0D0
  bf( 2) = + 1.0D0/36.0D0
  bf( 3) = - 1.0D0/3600.0D0
  bf( 4) = + 1.0D0/211680.0D0
  bf( 5) = - 1.0D0/10886400.0D0
  bf( 6) = + 1.0D0/526901760.0D0
  bf( 7) = - 4.064761645144226D-11
  bf( 8) = + 8.921691020456453D-13
  bf( 9) = - 1.993929586072108D-14
  bf(10) = + 4.518980029619918D-16

  rz = real(z)
  iz = aimag(z)
  nz = rz**2 + iz**2

  ! special cases
  if (iz .eq. 0) then
     if (rz .le. 1) cdli2 = dcmplx(dli2(rz), 0)
     if (rz .gt. 1) cdli2 = dcmplx(dli2(rz), -PI*log(rz))
     return
  elseif (nz .lt. EPSILON(1D0)) then
     cdli2 = z
     return
  endif

  ! transformation to |z| < 1, Re(z) <= 0.5
  if (rz .le. 0.5D0) then
     if (nz .gt. 1) then
        cy = -0.5D0*fast_cdlog(-z)**2 - PI*PI/6
        cz = -fast_cdlog(1 - 1/z)
        sgn = -1
     else ! nz <= 1
        cy = 0
        cz = -fast_cdlog(1 - z)
        sgn = 1
     endif
  else ! rz > 0.5D0
     if (nz .le. 2*rz) then
        cz = -fast_cdlog(z)
        cy = cz*fast_cdlog(1 - z) + PI*PI/6
        sgn = -1
     else ! nz > 2*rz
        cy = -0.5D0*fast_cdlog(-z)**2 - PI*PI/6
        cz = -fast_cdlog(1 - 1/z)
        sgn = -1
     endif
  endif

  cz2 = cz**2
  sum =                       &
      cz +                    &
      cz2 * (bf( 1) +         &
      cz  * (bf( 2) +         &
      cz2 * (bf( 3) +         &
      cz2 * (bf( 4) +         &
      cz2 * (bf( 5) +         &
      cz2 * (bf( 6) +         &
      cz2 * (bf( 7) +         &
      cz2 * (bf( 8) +         &
      cz2 * (bf( 9) +         &
      cz2 * (bf(10)))))))))))

  cdli2 = sgn*sum + cy

end function cdli2
