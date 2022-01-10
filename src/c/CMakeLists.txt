add_library(polylog_c
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

target_include_directories(polylog_c PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
