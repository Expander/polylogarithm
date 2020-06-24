subroutine dli2_wrapper(x, res) bind(C)
  use, intrinsic :: iso_c_binding
  implicit none
  real(c_double), intent(in)  :: x
  real(c_double), intent(out) :: res
  double precision dli2
  res = dli2(x)
end subroutine dli2_wrapper


subroutine cdli2_wrapper(re, im, res_re, res_im) bind(C)
  use, intrinsic :: iso_c_binding
  implicit none
  real(c_double), intent(in)  :: im, re
  real(c_double), intent(out) :: res_re, res_im
  double complex res, cdli2
  res = cdli2(dcmplx(re, im))
  res_re = real(res)
  res_im = aimag(res)
end subroutine cdli2_wrapper
