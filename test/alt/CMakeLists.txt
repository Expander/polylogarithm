add_subdirectory(algorithm_327)
add_subdirectory(algorithm_490)
add_subdirectory(babar)
add_subdirectory(bernoulli)
add_subdirectory(cephes)
add_subdirectory(hassani)
add_subdirectory(hdecay)
add_subdirectory(hollik)
add_subdirectory(koelbig)
add_subdirectory(feynhiggs)
add_subdirectory(morris)
add_subdirectory(pade)
add_subdirectory(sherpa)
add_subdirectory(spheno)
add_subdirectory(tsil)
add_subdirectory(wu)

add_library(alt alt.c)

target_include_directories(alt
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(alt
  PUBLIC
    algorithm_327
    algorithm_490
    babar
    bernoulli
    cephes
    hassani
    hdecay
    hollik
    koelbig
    feynhiggs
    morris
    pade
    sherpa
    spheno
    tsil
    wu
)
