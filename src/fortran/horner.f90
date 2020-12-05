!*********************************************************************
! This file is part of Polylogarithm.
!
! Polylogarithm is licenced under the MIT License.
!*********************************************************************


!*********************************************************************
!> @brief Evaluation of polynomial P(x) with len coefficients c
!> @param x real argument of P
!> @param c coefficients of P(x)
!> @param len number of coefficients
!> @return P(x)
!*********************************************************************

double precision function dhorner(x, c, len)
  implicit none
  integer :: len, i
  double precision :: x, c(len)

  dhorner = 0

  do i = len, 1, -1
     dhorner = dhorner*x + c(i)
  end do

end function dhorner
