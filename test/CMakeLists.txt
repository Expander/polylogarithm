find_package(GSL QUIET)

add_definitions(-DTEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")
add_subdirectory(algorithm_327)
add_subdirectory(algorithm_490)
add_subdirectory(cephes)
add_subdirectory(hollik)
add_subdirectory(tsil)

function(add_polylogarithm_test name)
  add_executable(
    ${name}
    ${name}.cpp
    )
  target_link_libraries(
    ${name}
    polylog
    ${ARGN}
    )
  add_test(NAME ${name} COMMAND ${name})
endfunction()

if(GSL_FOUND)
  add_polylogarithm_test(bench_Li GSL::gsl tsil algorithm_327 algorithm_490 cephes hollik)
endif()
  add_polylogarithm_test(test_Cl)
  add_polylogarithm_test(test_Cl1)
  add_polylogarithm_test(test_Cl2)
  add_polylogarithm_test(test_Cl3)
  add_polylogarithm_test(test_Cl4)
  add_polylogarithm_test(test_Cl5)
  add_polylogarithm_test(test_Cl6)
  add_polylogarithm_test(test_Li)
  add_polylogarithm_test(test_Li0)
  add_polylogarithm_test(test_Li1)
if(GSL_FOUND)
  add_polylogarithm_test(test_Li2 GSL::gsl tsil algorithm_327 algorithm_490 cephes hollik)
endif()
  add_polylogarithm_test(test_Li3 tsil)
  add_polylogarithm_test(test_Li4)
  add_polylogarithm_test(test_Li5)
  add_polylogarithm_test(test_Li6)
  add_polylogarithm_test(test_version)
