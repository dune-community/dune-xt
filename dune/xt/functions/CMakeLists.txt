# ~~~
# This file is part of the dune-xt-functions project:
#   https://github.com/dune-community/dune-xt-functions
# Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2016 - 2017)
#   Rene Milk       (2016, 2018)
#   Tobias Leibner  (2016)
#
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# ~~~

set(lib_dune_xt_functions_sources expression/mathexpr.cc)

if(DUNE_XT_WITH_PYTHON_BINDINGS)
  list(APPEND lib_dune_xt_functions_sources
       checkerboard.lib/alu_2d_simplex_conforming.cc
       checkerboard.lib/yasp_1d_equidistant_offset.cc
       checkerboard.lib/yasp_2d_equidistant_offset.cc
       checkerboard.lib/yasp_3d_equidistant_offset.cc
       constant.lib/alu_2d_simplex_conforming.cc
       constant.lib/yasp_1d_equidistant_offset.cc
       constant.lib/yasp_2d_equidistant_offset.cc
       constant.lib/yasp_3d_equidistant_offset.cc
       expression.lib/alu_2d_simplex_conforming.cc
       expression.lib/yasp_1d_equidistant_offset.cc
       expression.lib/yasp_2d_equidistant_offset.cc
       expression.lib/yasp_3d_equidistant_offset.cc
       interfaces.lib/alu_2d_simplex_conforming.cc
       interfaces.lib/yasp_1d_equidistant_offset.cc
       interfaces.lib/yasp_2d_equidistant_offset.cc
       interfaces.lib/yasp_3d_equidistant_offset.cc
       spe10/model1.lib/alu_2d_simplex_conforming.cc
       spe10/model1.lib/yasp_2d_equidistant_offset.cc
       spe10/model2.lib/yasp_3d_equidistant_offset.cc)
endif()
dune_library_add_sources(dunextfunctions SOURCES ${lib_dune_xt_functions_sources})

add_subdirectory(test EXCLUDE_FROM_ALL)

