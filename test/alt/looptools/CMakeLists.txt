if (ENABLE_FORTRAN)
  add_library(looptools
    ffcli2.F
    ffinit.F
    ffxli2.F
    ini.F
    Li2.F
    )
  add_compile_definitions(QUAD=0)
  target_include_directories(looptools PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
else()
  add_library(looptools INTERFACE)
endif()
