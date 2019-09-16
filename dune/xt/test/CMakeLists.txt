
enable_testing()

begin_testcases(dunextgrid dunextcommon dunextla)

end_testcases()

add_subdirectory(common)
add_subdirectory(grid)
add_subdirectory(la)
add_subdirectory(functions)

if (ALBERTA_FOUND)
  add_dune_alberta_flags(GRIDDIM 2 test_dd_glued_2d)
  add_dune_alberta_flags(GRIDDIM 3 test_dd_glued_3d)

  foreach(test ${dxt_test_binaries})
    if(${test} MATCHES alberta_1d)
      add_dune_alberta_flags(GRIDDIM 1 ${test})
    elseif(${test} MATCHES alberta_2d)
      add_dune_alberta_flags(GRIDDIM 2 ${test})
    elseif(${test} MATCHES alberta_3d)
      add_dune_alberta_flags(GRIDDIM 3 ${test})
    endif()
  endforeach()

  foreach(test ${dxt_test_binaries})
    if(${test} MATCHES 2d_simplex_alberta)
      add_dune_alberta_flags(GRIDDIM 2 ${test})
    elseif(${test} MATCHES 3d_simplex_alberta)
      add_dune_alberta_flags(GRIDDIM 3 ${test})
    endif()
  endforeach()

endif()
