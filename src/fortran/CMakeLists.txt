if(CMAKE_Fortran_COMPILER)
  add_library(polylogarithm_fortran
    fortran_wrappers.f90
    fast_clog.f90
    Cl2.f90
    Cl3.f90
    Cl4.f90
    Cl5.f90
    Cl6.f90
    Li2.f90
    Li3.f90
    Li4.f90
    Li5.f90
    Li6.f90
    )
  target_include_directories(polylogarithm_fortran PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
else()
  add_library(polylogarithm_fortran INTERFACE)
endif()
