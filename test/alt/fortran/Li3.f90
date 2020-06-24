!*********************************************************************
! This file is part of Polylogarithm.
!
! Polylogarithm is licenced under the GNU Lesser General Public
! License (GNU LGPL) version 3.
!*********************************************************************


!*********************************************************************
!> @brief Complex trilogarithm \f$\mathrm{Li}_3(z)\f$
!> @param z complex argument
!> @return \f$\mathrm{Li}_3(z)\f$
!> @author Alexander Voigt
!*********************************************************************
double complex function cdli3(z)
  implicit none
  double complex :: z, u, u2, u4, u8, c0, c1, lmz, rest, fast_pos_cdlog
  double precision :: rz, iz, nz, pz, lnz, arg
  double precision, parameter :: PI    = 3.1415926535897932D0
  double precision, parameter :: zeta2 = 1.6449340668482264D0
  double precision, parameter :: zeta3 = 1.2020569031595943D0
  double precision, parameter :: ln2   = 6.9314718055994531D-1 ! ln(2)
  double precision, parameter :: ln23  = 3.3302465198892948D-1 ! ln(2)^3
  double precision, parameter :: bf(18) = (/         &
      1.0D0                , -3.0D0/8.0D0          , &
      17.0D0/216.0D0       , -5.0D0/576.0D0        , &
      1.296296296296296D-04,  8.101851851851851D-05, &
     -3.419357160853759D-06, -1.328656462585034D-06, &
      8.660871756109851D-08,  2.526087595532039D-08, &
     -2.144694468364064D-09, -5.140110622012978D-10, &
      5.249582114600829D-11,  1.088775440663631D-11, &
     -1.277939609449369D-12, -2.369824177308745D-13, &
      3.104357887965462D-14,  5.261758629912506D-15  /)
  double precision, parameter :: cs(7) = (/          &
     -3.472222222222222D-03, 1.157407407407407D-05,  &
     -9.841899722852104D-08, 1.148221634332745D-09,  &
     -1.581572499080917D-11, 2.419500979252515D-13,  &
     -3.982897776989488D-15                          /)

  rz = real(z)
  iz = aimag(z)

   if (iz .eq. 0) then
      if (rz .eq. 0) then
         cdli3 = 0
         return
      endif
      if (rz .eq. 1) then
         cdli3 = zeta3
         return
      endif
      if (rz .eq. -1) then
         cdli3 = -0.75D0*zeta3
         return
      endif
      if (rz .eq. 0.5D0) then
         cdli3 = (-2*PI*PI*ln2 + 4*ln23 + 21*zeta3)/24
         return
      endif
   endif

   nz  = rz**2 + iz**2
   pz  = datan2(iz, rz)
   lnz = 0.5D0*log(nz)

   if (lnz*lnz + pz*pz .lt. 1) then ! |log(z)| < 1
      u = dcmplx(lnz, pz) ! log(z)
      u2 = u**2
      u4 = u2**2
      u8 = u4**2
      c0 = zeta3 + u*(zeta2 - u2/12)
      c1 = 0.25D0*(3 - 2*fast_pos_cdlog(-u))
      cdli3 =                                            &
         c0 +                                            &
         c1*u2 +                                         &
         u4*(cs(1) + u2*cs(2)) +                         &
         u8*(cs(3) + u2*cs(4) + u4*(cs(5) + u2*cs(6))) + &
         u8*u8*cs(7)
      return
   endif

   if (nz .le. 1) then
      u = -fast_pos_cdlog(1 - z)
      rest = 0
   else ! nz > 1
      if (pz .gt. 0) then
         arg = pz - PI
      else
         arg = pz + PI
      endif
      lmz = dcmplx(lnz, arg) ! log(-z)
      u = -fast_pos_cdlog(1 - 1/z)
      rest = -lmz*(lmz*lmz/6 + zeta2)
   endif

   u2 = u**2
   u4 = u2**2
   u8 = u4**2

   cdli3 =                                                   &
      rest +                                                 &
      u*bf(1) +                                              &
      u2*(bf(2) + u*bf(3)) +                                 &
      u4*(bf(4) + u*bf(5) + u2*(bf(6) + u*bf(7))) +          &
      u8*(bf(8) + u*bf(9) + u2*(bf(10) + u*bf(11)) +         &
          u4*(bf(12) + u*bf(13) + u2*(bf(14) + u*bf(15)))) + &
      u8*u8*(bf(16) + u*bf(17) + u2*bf(18))

end function cdli3
