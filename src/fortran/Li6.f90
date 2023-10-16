!*********************************************************************
! This file is part of Polylogarithm.
!
! Polylogarithm is licenced under the MIT License.
!*********************************************************************


!*********************************************************************
!> @brief Complex polylogarithm \f$\operatorname{Li}_6(z)\f$
!> @param z complex argument
!> @return \f$\operatorname{Li}_6(z)\f$
!> @author Alexander Voigt
!*********************************************************************
double complex function cdli6(z)
  implicit none
  double complex :: z, u, u2, u4, u8, c5, lmz, r, fast_pos_cdlog
  double precision :: rz, iz, nz, pz, lnz, arg, sgn
  double precision, parameter :: PI    = 3.1415926535897932D0
  double precision, parameter :: PI2   = 9.8696044010893586D0
  double precision, parameter :: PI4   = 97.409091034002437D0
  double precision, parameter :: PI6   = 961.38919357530444D0
  double precision, parameter :: zeta6 = 1.0173430619844491D0
  double precision, parameter :: c1    = 1.0369277551433699D0 ! zeta(5)
  double precision, parameter :: c2    = 0.54116161685556910D0
  double precision, parameter :: c3    = 0.20034281719326571D0
  double precision, parameter :: c4    = 0.068538919452009435D0
  double precision, parameter :: c6    = -1D0/1440
  double precision, parameter :: bf(18) = (/           &
      1.0D0                 , -31.0D0/64.0D0         , &
      1.5241340877914952D-01, -3.4365555877057613D-02, &
      5.7174797239368999D-03, -6.8180453746570645D-04, &
      4.9960361948734493D-05, -4.9166051196039048D-07, &
     -3.0632975161302164D-07,  1.4414599270849095D-08, &
      3.7272438230924107D-09, -3.7300867345487607D-10, &
     -5.1246526816085832D-11,  9.0541930956636683D-12, &
      6.7381882615512517D-13, -2.1215831150303135D-13, &
     -6.8408811719011698D-15,  4.8691178462005581D-15  /)
  double precision, parameter :: cs(5) = (/            &
     -1.6534391534391534D-05, 2.2964432686654909D-08,  &
     -9.9413128513657614D-11, 6.6912682653423394D-13,  &
     -5.7933058574392549D-15                           /)

  rz = real(z)
  iz = aimag(z)

   if (iz .eq. 0) then
      if (rz .eq. 0) then
         cdli6 = 0
         return
      endif
      if (rz .eq. 1) then
         cdli6 = zeta6
         return
      endif
      if (rz .eq. -1) then
         cdli6 = -31*zeta6/32
         return
      endif
   endif

   nz  = hypot(rz, iz)
   pz  = datan2(iz, rz)
   lnz = log(nz)

   if (lnz**2 + pz**2 .lt. 1) then ! |log(z)| < 1
      u = dcmplx(lnz, pz) ! log(z)
      u2 = u**2
      c5 = (137D0/60 - fast_pos_cdlog(-u))/120
      cdli6 = zeta6 + u * c1 + &
         u2 * (c2 + u * c3 +   &
         u2 * (c4 + u * c5 +   &
         u2 * (c6 +            &
         u * (cs(1) +          &
         u2 * (cs(2) +         &
         u2 * (cs(3) +         &
         u2 * (cs(4) +         &
         u2 * (cs(5)))))))))
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
      r = -31*PI6/15120 + lmz**2*(-7*PI4/720 + lmz**2*(-PI2/144 - lmz**2/720))
      sgn = -1
   endif

   u2 = u**2
   u4 = u2**2
   u8 = u4**2

   cdli6 =                                                      &
      r + sgn * (                                               &
         u*bf(1) +                                              &
         u2*(bf(2) + u*bf(3)) +                                 &
         u4*(bf(4) + u*bf(5) + u2*(bf(6) + u*bf(7))) +          &
         u8*(bf(8) + u*bf(9) + u2*(bf(10) + u*bf(11)) +         &
             u4*(bf(12) + u*bf(13) + u2*(bf(14) + u*bf(15)))) + &
         u8*u8*(bf(16) + u*bf(17) + u2*bf(18))                  &
      )

end function cdli6
