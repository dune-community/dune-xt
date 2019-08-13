# ~~~
# This file is part of the dune-xt-functions project:
#   https://github.com/dune-community/dune-xt-functions
# Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   René Fritze (2018)
#   Tim Keil    (2018)
# ~~~

dune_pybindxi_add_module(_function_interface_1d EXCLUDE_FROM_ALL function-interface-1d.cc)
dune_pybindxi_add_module(_function_interface_2d EXCLUDE_FROM_ALL function-interface-2d.cc)
dune_pybindxi_add_module(_function_interface_3d EXCLUDE_FROM_ALL function-interface-3d.cc)
dune_pybindxi_add_module(_gridfunction_interface_1d EXCLUDE_FROM_ALL gridfunction-interface-1d.cc)
dune_pybindxi_add_module(_gridfunction_interface_2d EXCLUDE_FROM_ALL gridfunction-interface-2d.cc)
dune_pybindxi_add_module(_gridfunction_interface_3d EXCLUDE_FROM_ALL gridfunction-interface-3d.cc)
dune_pybindxi_add_module(_checkerboard EXCLUDE_FROM_ALL checkerboard.cc)
dune_pybindxi_add_module(_constant EXCLUDE_FROM_ALL constant.cc)
dune_pybindxi_add_module(_expression EXCLUDE_FROM_ALL expression.cc)
dune_pybindxi_add_module(_function_as_grid_function EXCLUDE_FROM_ALL function-as-grid-function.cc)
dune_pybindxi_add_module(_indicator EXCLUDE_FROM_ALL indicator.cc)
dune_pybindxi_add_module(_spe10 EXCLUDE_FROM_ALL spe10.cc)
