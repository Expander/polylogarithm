if (ENABLE_FORTRAN)
  add_library(fortran
    cwrappers.f90
    fast_cdlog.f90
    horner.f90
    Li2.f90
    Li3.f90
    )

  target_include_directories(fortran PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
else()
  add_library(fortran INTERFACE)
endif()
