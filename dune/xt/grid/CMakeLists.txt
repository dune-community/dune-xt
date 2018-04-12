# ~~~
# This file is part of the dune-xt-grid project:
#   https://github.com/dune-community/dune-xt-grid
# Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2016 - 2017)
#   Rene Milk       (2016, 2018)
#   Tobias Leibner  (2016)
#
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# ~~~

set(lib_dune_xt_grid_sources boundaryinfo/interfaces.cc)

if(DUNE_XT_WITH_PYTHON_BINDINGS)
  list(APPEND lib_dune_xt_grid_sources
              gridprovider/cube.lib/alu_2d_simplex_conforming.cc
              gridprovider/cube.lib/yasp_1d_equidistant_offset.cc
              gridprovider/cube.lib/yasp_2d_equidistant_offset.cc
              gridprovider/cube.lib/yasp_3d_equidistant_offset.cc
              gridprovider/provider.lib/alu_2d_simplex_conforming.cc
              gridprovider/provider.lib/yasp_1d_equidistant_offset.cc
              gridprovider/provider.lib/yasp_2d_equidistant_offset.cc
              gridprovider/provider.lib/yasp_3d_equidistant_offset.cc
              walker.lib/alu_2d_simplex_conforming.cc
              walker.lib/yasp_1d_equidistant_offset.cc
              walker.lib/yasp_2d_equidistant_offset.cc
              walker.lib/yasp_3d_equidistant_offset.cc)
endif()
dune_library_add_sources(dunextgrid SOURCES ${lib_dune_xt_grid_sources})
add_dune_alberta_flags(GRIDDIM 1 SOURCE_ONLY dunextgrid)

include_directories(SYSTEM ${DUNE_XT_COMMON_TEST_DIR}/gtest)
add_subdirectory(test EXCLUDE_FROM_ALL)
