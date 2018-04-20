cmake_minimum_required(VERSION 3.1)
project(polylogarithm CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS OFF)

enable_testing()
add_subdirectory(src)
add_subdirectory(test)
