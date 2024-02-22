add_library(polylogarithm_c
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
set(C_HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/Cl2.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Cl3.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Cl4.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Cl5.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Cl6.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Li2.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Li3.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Li4.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Li5.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Li6.h
)
set_target_properties(polylogarithm_c PROPERTIES PUBLIC_HEADER "${C_HEADERS}")
