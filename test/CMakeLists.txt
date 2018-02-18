find_package(GSL REQUIRED)

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

add_dilog_test(test_dilog)
