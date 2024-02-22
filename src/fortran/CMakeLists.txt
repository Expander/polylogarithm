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
  set(FORTRAN_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/Cl2.fh
    ${CMAKE_CURRENT_SOURCE_DIR}/Cl3.fh
    ${CMAKE_CURRENT_SOURCE_DIR}/Cl4.fh
    ${CMAKE_CURRENT_SOURCE_DIR}/Cl5.fh
    ${CMAKE_CURRENT_SOURCE_DIR}/Cl6.fh
    ${CMAKE_CURRENT_SOURCE_DIR}/Li2.fh
    ${CMAKE_CURRENT_SOURCE_DIR}/Li3.fh
    ${CMAKE_CURRENT_SOURCE_DIR}/Li4.fh
    ${CMAKE_CURRENT_SOURCE_DIR}/Li5.fh
    ${CMAKE_CURRENT_SOURCE_DIR}/Li6.fh
    )
  set_target_properties(polylogarithm_c PROPERTIES PUBLIC_HEADER "${FORTRAN_HEADERS}")
else()
  add_library(polylogarithm_fortran INTERFACE)
endif()
