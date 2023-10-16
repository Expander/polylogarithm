!*********************************************************************
! This file is part of Polylogarithm.
!
! Polylogarithm is licenced under the MIT License.
!*********************************************************************


! Li_4(x) for x in [-1,0]
double precision function dli4_neg(x)
  implicit none
  double precision :: x, x2, x4, p, q
  double precision, parameter :: cp(6) = (/                &
      0.9999999999999999952D+0, -1.8532099956062184217D+0, &
      1.1937642574034898249D+0, -3.1817912243893560382D-1, &
      3.2268284189261624841D-2, -8.3773570305913850724D-4 /)
  double precision, parameter :: cq(7) = (/                &
      1.0000000000000000000D+0, -1.9157099956062165688D+0, &
      1.3011504531166486419D+0, -3.7975653506939627186D-1, &
      4.5822723996558783670D-2, -1.8023912938765272341D-3, &
      1.0199621542882314929D-5                            /)

  x2 = x*x
  x4 = x2*x2
  p = cp(1) + x*cp(2) + x2*(cp(3) + x*cp(4)) +             &
      x4*(cp(5) + x*cp(6))
  q = cq(1) + x*cq(2) + x2*(cq(3) + x*cq(4)) +             &
      x4*(cq(5) + x*cq(6) + x2*cq(7))

  dli4_neg = x*p/q

end function dli4_neg


! Li_4(x) for x in [0,1/2]
double precision function dli4_half(x)
  implicit none
  double precision :: x, x2, x4, p, q
  double precision, parameter :: cp(6) = (/                &
      1.0000000000000000414D+0, -2.0588072418045364525D+0, &
      1.4713328756794826579D+0, -4.2608608613069811474D-1, &
      4.2975084278851543150D-2, -6.8314031819918920802D-4 /)
  double precision, parameter :: cq(6) = (/                &
      1.0000000000000000000D+0, -2.1213072418045207223D+0, &
      1.5915688992789175941D+0, -5.0327641401677265813D-1, &
      6.1467217495127095177D-2, -1.9061294280193280330D-3 /)

  x2 = x*x
  x4 = x2*x2
  p = cp(1) + x*cp(2) + x2*(cp(3) + x*cp(4)) +             &
      x4*(cp(5) + x*cp(6))
  q = cq(1) + x*cq(2) + x2*(cq(3) + x*cq(4)) +             &
      x4*(cq(5) + x*cq(6))

  dli4_half = x*p/q

end function dli4_half


! Li_4(x) for x in [1/2,8/10]
double precision function dli4_mid(x)
  implicit none
  double precision :: x, x2, x4, p, q
  double precision, parameter :: cp(7) = (/                &
       3.2009826406098890447D-9, 9.9999994634837574160D-1, &
      -2.9144851228299341318D+0, 3.1891031447462342009D+0, &
      -1.6009125158511117090D+0, 3.5397747039432351193D-1, &
      -2.5230024124741454735D-2                           /)
  double precision, parameter :: cq(7) = (/                &
      1.0000000000000000000D+0, -2.9769855248411488460D+0, &
      3.3628208295110572579D+0, -1.7782471949702788393D+0, &
      4.3364007973198649921D-1, -3.9535592340362510549D-2, &
      5.7373431535336755591D-4                            /)

  x2 = x*x
  x4 = x2*x2
  p = cp(1) + x*cp(2) + x2*(cp(3) + x*cp(4)) +             &
      x4*(cp(5) + x*cp(6) + x2*cp(7))
  q = cq(1) + x*cq(2) + x2*(cq(3) + x*cq(4)) +             &
      x4*(cq(5) + x*cq(6) + x2*cq(7))

  dli4_mid = p/q

end function dli4_mid


! Li_4(x) for x in [8/10,1]
double precision function dli4_one(x)
  implicit none
  double precision :: x, l, l2
  double precision, parameter :: zeta2 = 1.6449340668482264D0
  double precision, parameter :: zeta3 = 1.2020569031595943D0
  double precision, parameter :: zeta4 = 1.0823232337111382D0

  l = log(x)
  l2 = l**2

  dli4_one = zeta4 + l*(zeta3 + l*(0.5D0*zeta2 + l*(11.0D0/36 &
     - 1.0D0/6*log(-l) + l*(-1.0D0/48 + l*(-1.0D0/1440        &
     + l2*(1.0D0/604800 - 1.0D0/91445760*l2))))))

end function dli4_one


!*********************************************************************
!> @brief Real 4-th order polylogarithm \f$\operatorname{Li}_4(x)\f$
!> @param x real argument
!> @return \f$\operatorname{Li}_4(x)\f$
!> @author Alexander Voigt
!*********************************************************************
double precision function dli4(x)
  implicit none
  double precision :: x, app, rest, sgn, l, l2
  double precision :: dli4_neg, dli4_half, dli4_mid, dli4_one
  double precision, parameter :: zeta2 = 1.6449340668482264D0
  double precision, parameter :: zeta4 = 1.0823232337111382D0

  ! transform x to [-1,1]
  if (x .lt. -1) then
     l = log(-x)
     l2 = l**2
     x = 1/x
     rest = -7.0D0/4*zeta4 + l2*(-0.5D0*zeta2 - 1.0D0/24*l2)
     sgn = -1
  elseif (x .eq. -1) then
     dli4 = -7.0D0/8*zeta4
     return
  elseif (x .lt. 1) then
     rest = 0
     sgn = 1
  elseif (x .eq. 1) then
     dli4 = zeta4
     return
  else ! x > 1
     l = log(x)
     l2 = l**2
     x = 1/x
     rest = 2*zeta4 + l2*(zeta2 - 1.0D0/24*l2)
     sgn = -1
  endif

  if (x .lt. 0) then
     app = dli4_neg(x)
  elseif (x .lt. 0.5D0) then
     app = dli4_half(x)
  elseif (x .lt. 0.8D0) then
     app = dli4_mid(x)
  else ! x <= 1
     app = dli4_one(x)
  endif

  dli4 = rest + sgn*app

end function dli4


!*********************************************************************
!> @brief Complex polylogarithm \f$\operatorname{Li}_4(z)\f$
!> @param z complex argument
!> @return \f$\operatorname{Li}_4(z)\f$
!> @author Alexander Voigt
!*********************************************************************
double complex function cdli4(z)
  implicit none
  double complex :: z, u, u2, u4, u8, c3, lmz, r, fast_pos_cdlog
  double precision :: rz, iz, nz, pz, lnz, arg, sgn, dli4
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
      if (rz .le. 1) then
         cdli4 = dli4(rz)
         return
      else
         lnz = log(rz)
         cdli4 = dcmplx(dli4(rz), -1D0/6D0*PI*lnz**3)
         return
      endif
   endif

   nz  = hypot(rz, iz)
   pz  = datan2(iz, rz)
   lnz = log(nz)

   if (lnz**2 + pz**2 .lt. 1) then ! |log(z)| < 1
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
