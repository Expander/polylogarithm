add_library(sherpa
  sherpa.c
)

target_include_directories(sherpa PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
