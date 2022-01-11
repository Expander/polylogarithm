cmake_minimum_required(VERSION 3.2)
project(polylogarithm CXX C)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED YES)
set(CMAKE_C_EXTENSIONS OFF)

include(CheckLanguage)
check_language(Fortran)

if(CMAKE_Fortran_COMPILER)
  enable_language(Fortran)
  set(CMAKE_Fortran_STANDARD 90)
  set(CMAKE_Fortran_STANDARD_REQUIRED YES)
  set(CMAKE_Fortran_EXTENSIONS OFF)
  add_compile_definitions(ENABLE_FORTRAN=1)
else()
  message(STATUS "No Fortran support")
endif()

enable_testing()
add_subdirectory(src)
add_subdirectory(test)
