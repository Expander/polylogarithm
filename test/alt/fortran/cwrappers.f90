subroutine dli2_wrapper(x, res) bind(C, name="dli2_wrapper")
  implicit none
  double precision, intent(in)  :: x
  double precision, intent(out) :: res
  double precision dli2
  res = dli2(x)
end subroutine dli2_wrapper


subroutine cdli2_wrapper(re, im, res_re, res_im) bind(C, name="cdli2_wrapper")
  implicit none
  double precision, intent(in)  :: im, re
  double precision, intent(out) :: res_re, res_im
  double complex res, cdli2
  res = cdli2(dcmplx(re, im))
  res_re = real(res)
  res_im = aimag(res)
end subroutine cdli2_wrapper
