# ~~~
# This file is part of the dune-xt project:
#   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
# Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2019 - 2020)
#   René Fritze     (2018 - 2019)
#   Tim Keil        (2018, 2020)
#   Tobias Leibner  (2019 - 2020)
# ~~~

dune_pybindxi_add_module(_grid_boundaryinfo_alldirichlet EXCLUDE_FROM_ALL boundaryinfo/alldirichlet.cc)
dune_pybindxi_add_module(_grid_boundaryinfo_allneumann EXCLUDE_FROM_ALL boundaryinfo/allneumann.cc)
dune_pybindxi_add_module(_grid_boundaryinfo_allreflecting EXCLUDE_FROM_ALL boundaryinfo/allreflecting.cc)
dune_pybindxi_add_module(_grid_boundaryinfo_interfaces EXCLUDE_FROM_ALL boundaryinfo/interfaces.cc)
dune_pybindxi_add_module(_grid_boundaryinfo_normalbased EXCLUDE_FROM_ALL boundaryinfo/normalbased.cc)
dune_pybindxi_add_module(_grid_boundaryinfo_functionbased EXCLUDE_FROM_ALL boundaryinfo/functionbased.cc)
dune_pybindxi_add_module(_grid_boundaryinfo_types EXCLUDE_FROM_ALL boundaryinfo/types.cc)
dune_pybindxi_add_module(_grid_element EXCLUDE_FROM_ALL element.cc)
dune_pybindxi_add_module(_grid_filters_base EXCLUDE_FROM_ALL filters/base.cc)
dune_pybindxi_add_module(_grid_filters_element EXCLUDE_FROM_ALL filters/element.cc)
dune_pybindxi_add_module(_grid_functors_boundary_detector EXCLUDE_FROM_ALL functors/boundary-detector.cc)
# dune_pybindxi_add_module(_grid_functors_bounding_box EXCLUDE_FROM_ALL functors/bounding-box.cc)
dune_pybindxi_add_module(_grid_functors_interfaces EXCLUDE_FROM_ALL functors/interfaces.cc)
# dune_pybindxi_add_module(_grid_functors_refinement EXCLUDE_FROM_ALL functors/refinement.cc)
dune_pybindxi_add_module(_grid_gridprovider_cube EXCLUDE_FROM_ALL gridprovider/cube.cc)
dune_pybindxi_add_module(_grid_gridprovider_gmsh EXCLUDE_FROM_ALL gridprovider/gmsh.cc)
dune_pybindxi_add_module(_grid_intersection EXCLUDE_FROM_ALL intersection.cc)
dune_pybindxi_add_module(_grid_gridprovider_provider EXCLUDE_FROM_ALL gridprovider/provider.cc)
dune_pybindxi_add_module(_grid_traits EXCLUDE_FROM_ALL traits.cc)
dune_pybindxi_add_module(_grid_walker EXCLUDE_FROM_ALL walker.cc)
dune_pybindxi_add_module(_grid_dd_glued_gridprovider_cube EXCLUDE_FROM_ALL dd_glued_gridprovider/cube.cc)
dune_pybindxi_add_module(_grid_dd_glued_gridprovider_provider EXCLUDE_FROM_ALL dd_glued_gridprovider/provider.cc)
