add_library(polylog_c
  c_wrappers.c
  Li2.c
  Li3.c
  Li4.c
  Li5.c
  Li6.c
)

target_include_directories(polylog_c PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
