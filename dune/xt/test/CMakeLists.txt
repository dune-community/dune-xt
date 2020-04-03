# ~~~
# This file is part of the dune-xt project:
#   https://github.com/dune-community/dune-xt
# Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2014 - 2016)
#   René Fritze     (2012 - 2016, 2019)
#   Tobias Leibner  (2015 - 2016, 2019 - 2020)
# ~~~

enable_testing()

dependencycheck(${xtcommon})

include_directories(SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/gtest)
dune_library_add_sources(gtest_dune_xt SOURCES common.cxx gtest/gtest-all.cxx)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/gtest
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/xt/test/
        FILES_MATCHING
        PATTERN "*")
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/main.hxx DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/xt/test/)

# this should ideally move back into the grid subdir, along with the files
dune_symlink_to_source_files(FILES
                             gmsh_2d_simplices.msh
                             dgf_1d_interval.dgf
                             dgf_2d_interval.dgf
                             dgf_3d_interval.dgf
                             dgf_4d_interval.dgf)

add_subdir_tests(common)
add_subdir_tests(grid)
add_subdir_tests(la)
add_subdir_tests(functions)

finalize_test_setup()
