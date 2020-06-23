subroutine dli2_wrapper(x, res) bind(C, name="dli2_wrapper")
  implicit none
  double precision, intent(in)  :: x
  double precision, intent(out) :: res
  double precision dli2

  res = dli2(x)
end subroutine dli2_wrapper
