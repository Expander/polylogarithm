add_subdirectory(c)
add_subdirectory(cpp)
add_subdirectory(fortran)

add_library(polylog INTERFACE)
target_include_directories(polylog
  INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/c
    ${CMAKE_CURRENT_SOURCE_DIR}/cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/fortran
)
target_link_libraries(polylog
  INTERFACE
    polylog_c
    polylog_cpp
    polylog_fortran
)
