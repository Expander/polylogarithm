add_library(fortran
  cwrappers.f90
  Li2.f90
)

target_include_directories(fortran PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
