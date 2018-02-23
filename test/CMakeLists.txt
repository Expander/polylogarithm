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

  add_polylogarithm_test(test_Cl2)
  add_polylogarithm_test(test_Li0)
  add_polylogarithm_test(test_Li1)
  add_polylogarithm_test(test_Li2)
  add_polylogarithm_test(test_Li3)
  add_polylogarithm_test(test_Li4)
  add_polylogarithm_test(test_Li5)
  add_polylogarithm_test(test_Li6)
endif()
