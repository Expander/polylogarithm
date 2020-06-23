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
      l = log(1 - x)
      y = x/(x - 1)
      r = -0.5D0*l*l
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
      l = log(x)
      y = 1/x
      r = PI*PI/3 - 0.5D0*l*l
      s = -1
   endif

  z = y - 0.25D0

  p = dhorner(z, cp, 8)
  q = dhorner(z, cq, 8)

  dli2 = r + s*y*p/q
end function dli2
