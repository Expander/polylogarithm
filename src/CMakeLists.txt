add_subdirectory(c)
add_subdirectory(cpp)
add_subdirectory(fortran)

add_library(polylogarithm INTERFACE)
target_include_directories(polylogarithm
  INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/c
    ${CMAKE_CURRENT_SOURCE_DIR}/cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/fortran
)
target_link_libraries(polylogarithm
  INTERFACE
    polylogarithm_c
    polylogarithm_cpp
    polylogarithm_fortran
)
set(HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/version.h)
set_target_properties(polylogarithm
  PROPERTIES PUBLIC_HEADER "${HEADERS}")
install(
  TARGETS
    polylogarithm
    polylogarithm_c
    polylogarithm_cpp
    polylogarithm_fortran
  PUBLIC_HEADER 
  DESTINATION
    ${CMAKE_INSTALL_INCLUDEDIR}/polylogarithm
)
