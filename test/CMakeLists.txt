find_package(GSL QUIET)

if(GSL_FOUND)
  set(LIBGSL GSL::gsl)
  add_definitions(-DENABLE_GSL=1)
endif()

add_definitions(-DTEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")

add_subdirectory(alt)

add_library(polylogarithm_test
    c_wrappers.c
    fortran_wrappers.f90
  )

function(add_polylogarithm_exec name)
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} alt polylogarithm polylogarithm_test ${ARGN})
endfunction()

function(add_polylogarithm_benc name)
  add_polylogarithm_exec(${test} ${name} ${ARGN})
endfunction()

function(add_polylogarithm_test name)
  add_polylogarithm_exec(${test} ${name} ${ARGN})
  add_test(NAME ${name} COMMAND ${name})
endfunction()

add_polylogarithm_benc(bench_Li ${LIBGSL})
add_polylogarithm_benc(bench_Cl ${LIBGSL})
add_polylogarithm_benc(bench_Sl)
add_polylogarithm_test(test_Cl)
add_polylogarithm_test(test_Cl1)
add_polylogarithm_test(test_Cl2 ${LIBGSL})
add_polylogarithm_test(test_Cl3)
add_polylogarithm_test(test_Cl4)
add_polylogarithm_test(test_Cl5)
add_polylogarithm_test(test_Cl6)
add_polylogarithm_test(test_eta)
add_polylogarithm_test(test_example)
add_polylogarithm_test(test_factorial)
add_polylogarithm_test(test_harmonic)
add_polylogarithm_test(test_Li)
add_polylogarithm_test(test_Li2 ${LIBGSL})
add_polylogarithm_test(test_Li3)
add_polylogarithm_test(test_Li4)
add_polylogarithm_test(test_Li5)
add_polylogarithm_test(test_Li6)
add_polylogarithm_test(test_log)
add_polylogarithm_test(test_Sl)
add_polylogarithm_test(test_zeta)
add_polylogarithm_test(test_version)
