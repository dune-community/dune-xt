# This file is part of the dune-xt-grid project:
#   https://github.com/dune-community/dune-xt-grid
# Copyright 2009-2017 dune-xt-grid developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2016 - 2017)
#   Rene Milk       (2016)
#   Tobias Leibner  (2016)

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
add_subdirectory(test EXCLUDE_FROM_ALL )

if(dune-pybindxi_FOUND)
  dune_pybindxi_add_helper_lib(_grid_bindings_helper
                               boundaryinfo.bindings.cc
                               )
  dune_pybindxi_add_module(_grid EXCLUDE_FROM_ALL bindings.cc)
  target_link_dune_default_libraries(_grid)
  target_link_libraries(_grid _grid_bindings_helper)
  add_dune_alberta_flags(GRIDDIM 2 _grid)
  if(DUNE_XT_WITH_PYTHON_BINDINGS)
    add_custom_target(bindings ALL DEPENDS _grid)
  else()
    add_custom_target(bindings DEPENDS _grid)
  endif()
endif()

