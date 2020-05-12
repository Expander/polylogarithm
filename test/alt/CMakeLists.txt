add_subdirectory(algorithm_327)
add_subdirectory(algorithm_490)
add_subdirectory(babar)
add_subdirectory(c)
add_subdirectory(cephes)
add_subdirectory(hassani)
add_subdirectory(hollik)
add_subdirectory(sherpa)
add_subdirectory(spheno)
add_subdirectory(tsil)

add_library(alt alt.c)

target_include_directories(alt
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(alt
  PRIVATE
    algorithm_327
    algorithm_490
    babar
    c
    cephes
    hassani
    hollik
    sherpa
    spheno
    tsil
)
