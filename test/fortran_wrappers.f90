!*********************************************************************
! This file is part of Polylogarithm.
!
! Polylogarithm is licenced under the MIT License.
!*********************************************************************


subroutine cl2_fortran(x, res) bind(C)
  use, intrinsic :: iso_c_binding
  implicit none
  real(c_double), intent(in)  :: x
  real(c_double), intent(out) :: res
  double precision dcl2
  res = dcl2(x)
end subroutine cl2_fortran


subroutine cl3_fortran(x, res) bind(C)
  use, intrinsic :: iso_c_binding
  implicit none
  real(c_double), intent(in)  :: x
  real(c_double), intent(out) :: res
  double precision dcl3
  res = dcl3(x)
end subroutine cl3_fortran


subroutine cl4_fortran(x, res) bind(C)
  use, intrinsic :: iso_c_binding
  implicit none
  real(c_double), intent(in)  :: x
  real(c_double), intent(out) :: res
  double precision dcl4
  res = dcl4(x)
end subroutine cl4_fortran


subroutine cl5_fortran(x, res) bind(C)
  use, intrinsic :: iso_c_binding
  implicit none
  real(c_double), intent(in)  :: x
  real(c_double), intent(out) :: res
  double precision dcl5
  res = dcl5(x)
end subroutine cl5_fortran


subroutine cl6_fortran(x, res) bind(C)
  use, intrinsic :: iso_c_binding
  implicit none
  real(c_double), intent(in)  :: x
  real(c_double), intent(out) :: res
  double precision dcl6
  res = dcl6(x)
end subroutine cl6_fortran


subroutine li2_fortran(x, res) bind(C)
  use, intrinsic :: iso_c_binding
  implicit none
  real(c_double), intent(in)  :: x
  real(c_double), intent(out) :: res
  double precision dli2
  res = dli2(x)
end subroutine li2_fortran


subroutine li3_fortran(x, res) bind(C)
  use, intrinsic :: iso_c_binding
  implicit none
  real(c_double), intent(in)  :: x
  real(c_double), intent(out) :: res
  double precision dli3
  res = dli3(x)
end subroutine li3_fortran


subroutine li4_fortran(x, res) bind(C)
  use, intrinsic :: iso_c_binding
  implicit none
  real(c_double), intent(in)  :: x
  real(c_double), intent(out) :: res
  double precision dli4
  res = dli4(x)
end subroutine li4_fortran


subroutine cli2_fortran(re, im, res_re, res_im) bind(C)
  use, intrinsic :: iso_c_binding
  implicit none
  real(c_double), intent(in)  :: im, re
  real(c_double), intent(out) :: res_re, res_im
  double complex res, cdli2
  res = cdli2(dcmplx(re, im))
  res_re = real(res)
  res_im = aimag(res)
end subroutine cli2_fortran


subroutine cli3_fortran(re, im, res_re, res_im) bind(C)
  use, intrinsic :: iso_c_binding
  implicit none
  real(c_double), intent(in)  :: im, re
  real(c_double), intent(out) :: res_re, res_im
  double complex res, cdli3
  res = cdli3(dcmplx(re, im))
  res_re = real(res)
  res_im = aimag(res)
end subroutine cli3_fortran


subroutine cli4_fortran(re, im, res_re, res_im) bind(C)
  use, intrinsic :: iso_c_binding
  implicit none
  real(c_double), intent(in)  :: im, re
  real(c_double), intent(out) :: res_re, res_im
  double complex res, cdli4
  res = cdli4(dcmplx(re, im))
  res_re = real(res)
  res_im = aimag(res)
end subroutine cli4_fortran


subroutine cli5_fortran(re, im, res_re, res_im) bind(C)
  use, intrinsic :: iso_c_binding
  implicit none
  real(c_double), intent(in)  :: im, re
  real(c_double), intent(out) :: res_re, res_im
  double complex res, cdli5
  res = cdli5(dcmplx(re, im))
  res_re = real(res)
  res_im = aimag(res)
end subroutine cli5_fortran


subroutine cli6_fortran(re, im, res_re, res_im) bind(C)
  use, intrinsic :: iso_c_binding
  implicit none
  real(c_double), intent(in)  :: im, re
  real(c_double), intent(out) :: res_re, res_im
  double complex res, cdli6
  res = cdli6(dcmplx(re, im))
  res_re = real(res)
  res_im = aimag(res)
end subroutine cli6_fortran
