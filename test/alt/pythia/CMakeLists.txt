add_library(pythia
  pythia.c
)

target_include_directories(pythia PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
