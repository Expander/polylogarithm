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
set(CXX_HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/Cl.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Cl1.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Cl2.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Cl3.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Cl4.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Cl5.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Cl6.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Li.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Li2.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Li3.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Li4.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Li5.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Li6.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Sl.hpp
)
set_target_properties(polylogarithm_cpp PROPERTIES PUBLIC_HEADER "${CXX_HEADERS}")
