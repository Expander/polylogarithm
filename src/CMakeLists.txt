add_subdirectory(cpp)

add_library(polylog INTERFACE)
target_include_directories(polylog INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/cpp)
target_link_libraries(polylog INTERFACE polylog_cpp)
