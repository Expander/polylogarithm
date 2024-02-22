add_library(polylogarithm_c
  c_wrappers.c
  Cl2.c
  Cl3.c
  Cl4.c
  Cl5.c
  Cl6.c
  Li2.c
  Li3.c
  Li4.c
  Li5.c
  Li6.c
)

target_include_directories(polylogarithm_c PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
SET(C_HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/fast_clog.h
)
set_target_properties(polylogarithm_c PROPERTIES PUBLIC_HEADER "${C_HEADERS}")