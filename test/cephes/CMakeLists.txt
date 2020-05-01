add_library(cephes
  cephes.c
)

target_include_directories(cephes PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
