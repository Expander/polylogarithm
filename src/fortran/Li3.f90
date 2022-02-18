!*********************************************************************
! This file is part of Polylogarithm.
!
! Polylogarithm is licenced under the MIT License.
!*********************************************************************


! Li_3(x) for x in [-1,0]
double precision function dli3_neg(x)
  implicit none
  double precision :: x, x2, x4, p, q
  double precision, parameter :: cp(6) = (/                &
      0.9999999999999999795D+0, -2.0281801754117129576D+0, &
      1.4364029887561718540D+0, -4.2240680435713030268D-1, &
      4.7296746450884096877D-2, -1.3453536579918419568D-3 /)
  double precision, parameter :: cq(7) = (/                &
      1.0000000000000000000D+0, -2.1531801754117049035D+0, &
      1.6685134736461140517D+0, -5.6684857464584544310D-1, &
      8.1999463370623961084D-2, -4.0756048502924149389D-3, &
      3.4316398489103212699D-5                            /)

  x2 = x*x
  x4 = x2*x2
  p = cp(1) + x*cp(2) + x2*(cp(3) + x*cp(4)) +             &
       x4*(cp(5) + x*cp(6))
  q = cq(1) + x*cq(2) + x2*(cq(3) + x*cq(4)) +             &
       x4*(cq(5) + x*cq(6) + x2*cq(7))

  dli3_neg = x*p/q

end function dli3_neg


! Li_3(x) for x in [0,1/2]
double precision function dli3_pos(x)
  implicit none
  double precision :: x, x2, x4, p, q
  double precision, parameter :: cp(6) = (/                &
      0.9999999999999999893D+0, -2.5224717303769789628D+0, &
      2.3204919140887894133D+0, -9.3980973288965037869D-1, &
      1.5728950200990509052D-1, -7.5485193983677071129D-3 /)
  double precision, parameter :: cq(7) = (/                &
      1.0000000000000000000D+0, -2.6474717303769836244D+0, &
      2.6143888433492184741D+0, -1.1841788297857667038D+0, &
      2.4184938524793651120D-1, -1.8220900115898156346D-2, &
      2.4927971540017376759D-4                            /)

  x2 = x*x
  x4 = x2*x2
  p = cp(1) + x*cp(2) + x2*(cp(3) + x*cp(4)) +             &
       x4*(cp(5) + x*cp(6))
  q = cq(1) + x*cq(2) + x2*(cq(3) + x*cq(4)) +             &
       x4*(cq(5) + x*cq(6) + x2*cq(7))

  dli3_pos = x*p/q

end function dli3_pos


!*********************************************************************
!> @brief Real trilogarithm \f$\operatorname{Li}_3(x)\f$
!> @param x real argument
!> @return \f$\operatorname{Li}_3(x)\f$
!> @author Alexander Voigt
!*********************************************************************
double precision function dli3(x)
  implicit none
  double precision :: x, neg, pos, rest, sgn, l
  double precision :: dli3_neg, dli3_pos
  double precision, parameter :: zeta2 = 1.6449340668482264D0
  double precision, parameter :: zeta3 = 1.2020569031595943D0

  ! transformation to [-1,0] and [0,1/2]
  if (x .lt. -1) then
     l = log(-x)
     neg = dli3_neg(1/x)
     pos = 0
     sgn = 1
     rest = -l*(zeta2 + 1.0D0/6*l**2)
  elseif (x .eq. -1) then
     dli3 = -0.75D0*zeta3
     return
  elseif (x .lt. 0) then
     neg = dli3_neg(x)
     pos = 0
     sgn = 1
     rest = 0
  elseif (x .eq. 0) then
     dli3 = 0
     return
  elseif (x .lt. 0.5D0) then
     neg = 0
     pos = dli3_pos(x)
     sgn = 1
     rest = 0
  elseif (x .eq. 0.5D0) then
     dli3 = 0.53721319360804020D0
     return
  elseif (x .lt. 1) then
     l = log(x)
     neg = dli3_neg((x - 1)/x)
     pos = dli3_pos(1 - x)
     sgn = -1
     rest = zeta3 + l*(zeta2 + l*(-0.5D0*log(1 - x) + 1.0D0/6*l))
  elseif (x .eq. 1) then
     dli3 = zeta3
     return
  elseif (x .lt. 2) then
     l = log(x)
     neg = dli3_neg(1 - x)
     pos = dli3_pos((x - 1)/x)
     sgn = -1
     rest = zeta3 + l*(zeta2 + l*(-0.5D0*log(x - 1) + 1.0D0/6*l))
  else ! x >= 2.0D0
     l = log(x)
     neg = 0
     pos = dli3_pos(1/x)
     sgn = 1
     rest = l*(2*zeta2 - 1.0D0/6*l**2)
  endif

  dli3 = rest + sgn*(neg + pos)

end function dli3


!*********************************************************************
!> @brief Complex trilogarithm \f$\operatorname{Li}_3(z)\f$
!> @param z complex argument
!> @return \f$\operatorname{Li}_3(z)\f$
!> @author Alexander Voigt
!*********************************************************************
double complex function cdli3(z)
  implicit none
  double complex :: z, u, u2, u4, u8, c0, c1, lmz, rest, fast_pos_cdlog
  double precision :: rz, iz, nz, pz, lnz, arg
  double precision, parameter :: PI    = 3.1415926535897932D0
  double precision, parameter :: zeta2 = 1.6449340668482264D0
  double precision, parameter :: zeta3 = 1.2020569031595943D0
  double precision, parameter :: bf(18) = (/           &
      1.0D0                 , -3.0D0/8.0D0           , &
      17.0D0/216.0D0        , -5.0D0/576.0D0         , &
      1.2962962962962963D-04,  8.1018518518518519D-05, &
     -3.4193571608537595D-06, -1.3286564625850340D-06, &
      8.6608717561098513D-08,  2.5260875955320400D-08, &
     -2.1446944683640648D-09, -5.1401106220129789D-10, &
      5.2495821146008294D-11,  1.0887754406636318D-11, &
     -1.2779396094493695D-12, -2.3698241773087452D-13, &
      3.1043578879654623D-14,  5.2617586299125061D-15  /)
  double precision, parameter :: cs(7) = (/            &
     -3.4722222222222222D-03,  1.1574074074074074D-05, &
     -9.8418997228521038D-08,  1.1482216343327454D-09, &
     -1.5815724990809166D-11,  2.4195009792525152D-13, &
     -3.9828977769894877D-15                           /)

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
         cdli3 = 0.53721319360804020D0
         return
      endif
   endif

   nz  = rz**2 + iz**2
   pz  = datan2(iz, rz)
   lnz = 0.5D0*log(nz)

   if (lnz**2 + pz**2 .lt. 1) then ! |log(z)| < 1
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
      rest = -lmz*(lmz**2/6 + zeta2)
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
