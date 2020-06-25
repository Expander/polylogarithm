!*********************************************************************
! This file is part of Polylogarithm.
!
! Polylogarithm is licenced under the GNU Lesser General Public
! License (GNU LGPL) version 3.
!*********************************************************************


!*********************************************************************
!> @brief Complex polylogarithm \f$\mathrm{Li}_4(z)\f$
!> @param z complex argument
!> @return \f$\mathrm{Li}_4(z)\f$
!> @author Alexander Voigt
!*********************************************************************
double complex function cdli4(z)
  implicit none
  double complex :: z, u, u2, u4, u8, c3, lmz, r, fast_pos_cdlog
  double precision :: rz, iz, nz, pz, lnz, arg, sgn
  double precision, parameter :: PI    = 3.1415926535897932D0
  double precision, parameter :: PI2   = 9.8696044010893586D0
  double precision, parameter :: PI4   = 97.409091034002437D0
  double precision, parameter :: zeta4 = 1.0823232337111382D0
  double precision, parameter :: c1    = 1.2020569031595943D0
  double precision, parameter :: c2    = 0.82246703342411322D0
  double precision, parameter :: c4    = -1D0/48D0
  double precision, parameter :: bf(18) = (/           &
      1.0D0                 , -7.0D0/16.0D0          , &
      1.1651234567901235D-01, -1.9820601851851852D-02, &
      1.9279320987654321D-03, -3.1057098765432099D-05, &
     -1.5624009114857835D-05,  8.4851235467732066D-07, &
      2.2909616603189711D-07, -2.1832614218526917D-08, &
     -3.8828248791720156D-09,  5.4462921032203321D-10, &
      6.9608052106827254D-11, -1.3375737686445215D-11, &
     -1.2784852685266572D-12,  3.2605628580248922D-13, &
      2.3647571168618257D-14, -7.9231351220311617D-15  /)
  double precision, parameter :: cs(7) = (/            &
     -6.9444444444444444D-04, 1.6534391534391534D-06,  &
     -1.0935444136502338D-08, 1.0438378493934049D-10,  &
     -1.2165942300622435D-12, 1.6130006528350101D-14,  &
     -2.3428810452879340D-16                           /)

  rz = real(z)
  iz = aimag(z)

   if (iz .eq. 0) then
      if (rz .eq. 0) then
         cdli4 = 0
         return
      endif
      if (rz .eq. 1) then
         cdli4 = zeta4
         return
      endif
      if (rz .eq. -1) then
         cdli4 = -7*PI4/720
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
      c3 = (11D0/6 - fast_pos_cdlog(-u))/6
      cdli4 = zeta4 + u2 * (c2 + u2 * c4) +                   &
          u * (                                               &
              c1 +                                            &
              c3*u2 +                                         &
              u4*(cs(1) + u2*cs(2)) +                         &
              u8*(cs(3) + u2*cs(4) + u4*(cs(5) + u2*cs(6))) + &
              u8*u8*cs(7)                                     &
          )
      return
   endif

   if (nz .le. 1) then
      u = -fast_pos_cdlog(1 - z)
      r = 0
      sgn = 1
   else ! nz > 1
      if (pz .gt. 0) then
         arg = pz - PI
      else
         arg = pz + PI
      endif
      lmz = dcmplx(lnz, arg) ! log(-z)
      u = -fast_pos_cdlog(1 - 1/z)
      r = (-7*PI4 + lmz**2*(-30*PI2 - 15*lmz**2))/360
      sgn = -1
   endif

   u2 = u**2
   u4 = u2**2
   u8 = u4**2

   cdli4 =                                                      &
      r + sgn * (                                               &
         u*bf(1) +                                              &
         u2*(bf(2) + u*bf(3)) +                                 &
         u4*(bf(4) + u*bf(5) + u2*(bf(6) + u*bf(7))) +          &
         u8*(bf(8) + u*bf(9) + u2*(bf(10) + u*bf(11)) +         &
             u4*(bf(12) + u*bf(13) + u2*(bf(14) + u*bf(15)))) + &
         u8*u8*(bf(16) + u*bf(17) + u2*bf(18))                  &
      )

end function cdli4
