add_library(polylogarithm_cpp
  Cl.cpp
  Cl1.cpp
  Cl2.cpp
  Cl3.cpp
  Cl4.cpp
  Cl5.cpp
  Cl6.cpp
  eta.cpp
  factorial.cpp
  harmonic.cpp
  Li.cpp
  Li2.cpp
  Li3.cpp
  Li4.cpp
  Li5.cpp
  Li6.cpp
  Sl.cpp
  zeta.cpp
  )
target_include_directories(polylogarithm_cpp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
FILE(GLOB HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)
set_target_properties(polylogarithm_cpp PROPERTIES PUBLIC_HEADER "${HEADERS}")
