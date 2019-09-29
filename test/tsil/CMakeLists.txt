add_library(tsil
  dilog.c
  trilog.c
)

target_include_directories(tsil PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
