!*********************************************************************
! This file is part of Polylogarithm.
!
! Polylogarithm is licenced under the MIT License.
!*********************************************************************


!*********************************************************************
!> @brief Complex polylogarithm \f$\operatorname{Li}_5(z)\f$
!> @param z complex argument
!> @return \f$\operatorname{Li}_5(z)\f$
!> @author Alexander Voigt
!*********************************************************************
double complex function cdli5(z)
  implicit none
  double complex :: z, u, u2, u4, u8, c4, lmz, rest, pos_cdlog
  double precision :: rz, iz, nz, pz, lnz, arg
  double precision, parameter :: PI    = 3.1415926535897932D0
  double precision, parameter :: PI2   = 9.8696044010893586D0
  double precision, parameter :: PI4   = 97.409091034002437D0
  double precision, parameter :: zeta5 = 1.0369277551433699D0
  double precision, parameter :: c1    = 1.0823232337111382D0 ! zeta(4)
  double precision, parameter :: c2    = 0.60102845157979714D0 ! zeta(3)/2
  double precision, parameter :: c3    = 0.27415567780803774D0
  double precision, parameter :: c5    = -1D0/240
  double precision, parameter :: bf(19) = (/           &
      1.0D0                 , -15.0D0/32.0D0         , &
      1.3953189300411523D-01, -2.8633777006172840D-02, &
      4.0317412551440329D-03, -3.3985018004115226D-04, &
      4.5445184621617666D-06,  2.3916808048569012D-06, &
     -1.2762692600122747D-07, -3.1628984306505932D-08, &
      3.2848118445335192D-09,  4.7613713995660579D-10, &
     -8.0846898171909830D-11, -7.2387648587737207D-12, &
      1.9439760115173968D-12,  1.0256978405977236D-13, &
     -4.6180551009884830D-14, -1.1535857196470580D-15, &
      1.0903545401333394D-15                           /)
  double precision, parameter :: cs(6) = (/            &
     -1.1574074074074074D-04, 2.0667989417989418D-07,  &
     -1.0935444136502338D-09, 8.6986487449450412D-12,  &
     -8.6899587861588824D-14, 1.0081254080218813D-15   /)

  rz = real(z)
  iz = aimag(z)

   if (iz .eq. 0) then
      if (rz .eq. 0) then
         cdli5 = dcmplx(rz, iz)
         return
      endif
      if (rz .eq. 1) then
         cdli5 = dcmplx(zeta5, iz)
         return
      endif
      if (rz .eq. -1) then
         cdli5 = dcmplx(-15*zeta5/16, iz)
         return
      endif
   endif

   nz  = hypot(rz, iz)
   pz  = datan2(iz, rz)
   lnz = log(nz)

   if (lnz**2 + pz**2 .lt. 1) then ! |log(z)| < 1
      u = dcmplx(lnz, pz) ! log(z)
      u2 = u**2
      c4 = (25D0/12 - pos_cdlog(-u))/24
      cdli5 =                &
         zeta5 + u * c1 +    &
         u2 * (c2 + u * c3 + &
         u2 * (c4 + u * c5 + &
         u2 * (cs(1) +       &
         u2 * (cs(2) +       &
         u2 * (cs(3) +       &
         u2 * (cs(4) +       &
         u2 * (cs(5) +       &
         u2 * (cs(6)))))))))
      return
   endif

   if (nz .le. 1) then
      u = -pos_cdlog(1 - z)
      rest = 0
   else ! nz > 1
      if (pz .gt. 0) then
         arg = pz - PI
      else
         arg = pz + PI
      endif
      lmz = dcmplx(lnz, arg) ! log(-z)
      u = -pos_cdlog(1 - 1/z)
      rest = -lmz*(7*PI4 + lmz**2*(10*PI2 + 3*lmz**2))/360
   endif

   u2 = u**2
   u4 = u2**2
   u8 = u4**2

   cdli5 =                                                   &
      rest +                                                 &
      u*bf(1) +                                              &
      u2*(bf(2) + u*bf(3)) +                                 &
      u4*(bf(4) + u*bf(5) + u2*(bf(6) + u*bf(7))) +          &
      u8*(bf(8) + u*bf(9) + u2*(bf(10) + u*bf(11)) +         &
          u4*(bf(12) + u*bf(13) + u2*(bf(14) + u*bf(15)))) + &
      u8*u8*(bf(16) + u*bf(17) + u2*(bf(18) + u*bf(19)))

end function cdli5
