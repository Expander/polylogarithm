add_library(bernoulli
  bernoulli.c
)

target_include_directories(bernoulli PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
