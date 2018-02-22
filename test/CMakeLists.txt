find_package(GSL)

if(GSL_FOUND)
  function(add_dilog_test name)
    add_executable(
      ${name}
      ${name}.cpp
      )
    target_link_libraries(
      ${name}
      dilog
      GSL::gsl
      )
    add_test(NAME ${name} COMMAND ${name})
  endfunction()

  add_dilog_test(test_Cl2)
  add_dilog_test(test_Li2)
  add_dilog_test(test_Li3)
  add_dilog_test(test_Li4)
endif()
