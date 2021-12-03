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
!> Implemented as an economized Pade approximation with a
!> maximum error of 4.16e-18.
!*********************************************************************
      double precision function dli2(x)
      implicit none
      double precision x, y, r, s, z, z2, z4, p, q, l, PI
      parameter (PI = 3.14159265358979324D0)
      double precision cp(8), cq(8)
      data cp / 1.0706105563309304277D+0,
     &         -4.5353562730201404017D+0,
     &          7.4819657596286408905D+0,
     &         -6.0516124315132409155D+0,
     &          2.4733515209909815443D+0,
     &         -4.6937565143754629578D-1,
     &          3.1608910440687221695D-2,
     &         -2.4630612614645039828D-4 /
      data cq / 1.0000000000000000000D+0,
     &         -4.5355682121856044935D+0,
     &          8.1790029773247428573D+0,
     &         -7.4634190853767468810D+0,
     &          3.6245392503925290187D+0,
     &         -8.9936784740041174897D-1,
     &          9.8554565816757007266D-2,
     &         -3.2116618742475189569D-3 /

      ! transform to [0, 1/2)
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
          dli2 = 0
          return
       elseif (x .lt. 0.5D0) then
          y = x
          r = 0
          s = 1
       elseif (x .lt. 1) then
          y = 1 - x
          r = PI**2/6 - log(x)*log(1 - x)
          s = -1
       elseif (x .eq. 1) then
          dli2 = PI**2/6
          return
       elseif (x .lt. 2) then
          l = log(x)
          y = 1 - 1/x
          r = PI**2/6 - l*(log(1 - 1/x) + 0.5D0*l)
          s = 1
       else
          y = 1/x
          r = PI**2/3 - 0.5D0*log(x)**2
          s = -1
       endif

      z = y - 0.25D0
      z2 = z*z
      z4 = z2*z2

      p = cp(1) + z * cp(2) + z2 * (cp(3) + z * cp(4)) +
     &    z4 * (cp(5) + z * cp(6) + z2 * (cp(7) + z * cp(8)))
      q = cq(1) + z * cq(2) + z2 * (cq(3) + z * cq(4)) +
     &    z4 * (cq(5) + z * cq(6) + z2 * (cq(7) + z * cq(8)))

      dli2 = r + s*y*p/q

      end


!*********************************************************************
!> @brief Complex dilogarithm \f$\operatorname{Li}_2(z)\f$
!> @param z complex argument
!> @return \f$\operatorname{Li}_2(z)\f$
!> @note Implementation adapted from SPheno by Alexander Voigt
!*********************************************************************
      double complex function cdli2(z)
      implicit none
      double complex z, rest, u, u2, u4, sum, fast_cdlog
      double precision rz, iz, nz, sgn, dli2, PI
      parameter (PI = 3.14159265358979324D0)
      double precision bf(10)
      data bf /- 2.5000000000000000D-01,
     &         + 2.7777777777777778D-02,
     &         - 2.7777777777777778D-04,
     &         + 4.7241118669690098D-06,
     &         - 9.1857730746619636D-08,
     &         + 1.8978869988970999D-09,
     &         - 4.0647616451442255D-11,
     &         + 8.9216910204564526D-13,
     &         - 1.9939295860721076D-14,
     &         + 4.5189800296199182D-16 /

      rz = real(z)
      iz = dimag(z)

      ! special cases
      if (iz .eq. 0) then
         if (rz .le. 1) cdli2 = dcmplx(dli2(rz), 0)
         if (rz .gt. 1) cdli2 = dcmplx(dli2(rz), -PI*log(rz))
         return
      endif

      nz = rz**2 + iz**2

      if (nz .lt. 1.0D-15) then
         cdli2 = z*(1 + 0.25D0*z)
         return
      endif

      ! transformation to |z| < 1, Re(z) <= 0.5
      if (rz .le. 0.5D0) then
         if (nz .gt. 1) then
            u = -fast_cdlog(1 - 1/z)
            rest = -0.5D0*fast_cdlog(-z)**2 - PI**2/6
            sgn = -1
         else ! nz <= 1
            u = -fast_cdlog(1 - z)
            rest = 0
            sgn = 1
         endif
      else ! rz > 0.5D0
         if (nz .le. 2*rz) then
            u = -fast_cdlog(z)
            rest = u*fast_cdlog(1 - z) + PI**2/6
            sgn = -1
         else ! nz > 2*rz
            u = -fast_cdlog(1 - 1/z)
            rest = -0.5D0*fast_cdlog(-z)**2 - PI**2/6
            sgn = -1
         endif
      endif

      u2 = u**2
      u4 = u2**2
      sum =
     &   u +
     &   u2 * (bf(1) +
     &   u  * (bf(2) +
     &   u2 * (
     &       bf(3) +
     &       u2*bf(4) +
     &       u4*(bf(5) + u2*bf(6)) +
     &       u4*u4*(bf(7) + u2*bf(8) + u4*(bf(9) + u2*bf(10)))
     &   )))

      cdli2 = sgn*sum + rest

      end
