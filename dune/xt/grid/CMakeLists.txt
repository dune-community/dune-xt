# ~~~
# This file is part of the dune-xt project:
#   https://github.com/dune-community/dune-xt
# Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2013 - 2017)
#   René Fritze     (2012 - 2016, 2018 - 2019)
#   Tobias Leibner  (2016, 2019 - 2020)
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
dune_library_add_sources(dunext SOURCES ${lib_dune_xt_grid_sources})
