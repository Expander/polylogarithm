find_package(GSL QUIET)

if(GSL_FOUND)
  function(add_polylogarithm_test name)
    add_executable(
      ${name}
      ${name}.cpp
      )
    target_link_libraries(
      ${name}
      polylog
      GSL::gsl
      )
    add_test(NAME ${name} COMMAND ${name})
  endfunction()

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
  add_polylogarithm_test(test_Li2)
  add_polylogarithm_test(test_Li3)
  add_polylogarithm_test(test_Li4)
  add_polylogarithm_test(test_Li5)
  add_polylogarithm_test(test_Li6)
  add_polylogarithm_test(test_version)

  add_polylogarithm_test(bench_Li2)
  add_polylogarithm_test(bench_Li3)
endif()
