# ~~~
# This file is part of the dune-xt project:
#   https://github.com/dune-community/dune-xt
# Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2019)
#   René Fritze     (2018 - 2019)
#   Tim Keil        (2018)
#   Tobias Leibner  (2019 - 2020)
# ~~~

dune_pybindxi_add_module(_functions_checkerboard EXCLUDE_FROM_ALL checkerboard.cc)
dune_pybindxi_add_module(_functions_constant EXCLUDE_FROM_ALL constant.cc)
dune_pybindxi_add_module(_functions_expression EXCLUDE_FROM_ALL expression.cc)
dune_pybindxi_add_module(_functions_function_as_grid_function EXCLUDE_FROM_ALL function-as-grid-function.cc)
dune_pybindxi_add_module(_functions_function_interface_1d EXCLUDE_FROM_ALL function-interface-1d.cc)
dune_pybindxi_add_module(_functions_function_interface_2d EXCLUDE_FROM_ALL function-interface-2d.cc)
dune_pybindxi_add_module(_functions_function_interface_3d EXCLUDE_FROM_ALL function-interface-3d.cc)
dune_pybindxi_add_module(_functions_gridfunction EXCLUDE_FROM_ALL gridfunction.cc)
dune_pybindxi_add_module(_functions_indicator EXCLUDE_FROM_ALL indicator.cc)
dune_pybindxi_add_module(_functions_interfaces_grid_function_1d EXCLUDE_FROM_ALL interfaces/grid-function_1d.cc)
dune_pybindxi_add_module(_functions_interfaces_grid_function_2d EXCLUDE_FROM_ALL interfaces/grid-function_2d.cc)
dune_pybindxi_add_module(_functions_interfaces_grid_function_3d EXCLUDE_FROM_ALL interfaces/grid-function_3d.cc)
dune_pybindxi_add_module(_functions_inverse EXCLUDE_FROM_ALL inverse.cc)
dune_pybindxi_add_module(_functions_parametric_expression EXCLUDE_FROM_ALL parametric-expression.cc)
dune_pybindxi_add_module(_functions_spe10 EXCLUDE_FROM_ALL spe10.cc)
