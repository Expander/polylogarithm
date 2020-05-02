add_subdirectory(algorithm_327)
add_subdirectory(algorithm_490)
add_subdirectory(cephes)
add_subdirectory(hollik)
add_subdirectory(tsil)

add_library(alt alt.c)

target_include_directories(alt
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(alt
  PRIVATE
    algorithm_327
    algorithm_490
    cephes
    hollik
    tsil
)
