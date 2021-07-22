find_package(GSL QUIET)

if(GSL_FOUND)
  set(LIBGSL GSL::gsl)
  add_definitions(-DENABLE_GSL=1)
endif()

add_definitions(-DTEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")

add_subdirectory(alt)

function(add_polylogarithm_test name)
  add_executable(
    ${name}
    ${name}.cpp
    )
  target_link_libraries(
    ${name}
    alt
    polylog
    ${ARGN}
    )
  add_test(NAME ${name} COMMAND ${name})
endfunction()

add_polylogarithm_test(bench_Li ${LIBGSL})
add_polylogarithm_test(bench_Cl ${LIBGSL})
add_polylogarithm_test(test_Cl)
add_polylogarithm_test(test_Cl2 ${LIBGSL})
add_polylogarithm_test(test_Cl3)
add_polylogarithm_test(test_Cl4)
add_polylogarithm_test(test_Cl5)
add_polylogarithm_test(test_Cl6)
add_polylogarithm_test(test_Li)
add_polylogarithm_test(test_Li2 ${LIBGSL})
add_polylogarithm_test(test_Li3)
add_polylogarithm_test(test_Li4)
add_polylogarithm_test(test_Li5)
add_polylogarithm_test(test_Li6)
add_polylogarithm_test(test_version)
