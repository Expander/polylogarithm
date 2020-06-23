subroutine li2_fortran(x, res) bind(C, name="li2_fortran_")
  implicit none
  double precision, intent(in)  :: x
  double precision, intent(out) :: res
  double precision dli2

  res = dli2(x)
end subroutine li2_fortran
