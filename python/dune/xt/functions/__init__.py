# ~~~
# This file is part of the dune-xt project:
#   https://github.com/dune-community/dune-xt
# Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2017, 2019)
#   René Fritze     (2018 - 2019)
#   Tim Keil        (2018)
#   Tobias Leibner  (2019 - 2020)
# ~~~

from tempfile import NamedTemporaryFile

from dune.xt import guarded_import
from dune.xt.common.vtk.plot import plot

for mod_name in (
        '_functions_checkerboard',
        '_functions_constant',
        '_functions_divergence',
        '_functions_elementwise_diameter',
        '_functions_elementwise_minimum',
        '_functions_expression',
        '_functions_function_as_grid_function',
        '_functions_function_interface_1d',
        '_functions_function_interface_2d',
        '_functions_function_interface_3d',
        '_functions_gradient',
        '_functions_gridfunction',
        '_functions_indicator',
        '_functions_interfaces_element_function_1d',
        '_functions_interfaces_element_function_2d',
        '_functions_interfaces_element_function_3d',
        '_functions_interfaces_grid_function_1d',
        '_functions_interfaces_grid_function_2d',
        '_functions_interfaces_grid_function_3d',
        '_functions_inverse',
        '_functions_parametric_expression',
        '_functions_spe10',
):
    guarded_import(globals(), 'dune.xt.functions', mod_name)

from dune.xt.functions._functions_gridfunction import GridFunction


def visualize_function(function, grid, subsampling=False):
    assert function.dim_domain == 2, f'Not implemented yet for {function.dim_domain}-dimensional grids!'
    assert function.dim_range == 1, f'Not implemented yet for {function.dim_range}-dimensional functions!'
    tmpfile = NamedTemporaryFile(mode='wb', delete=False, suffix='.vtu').name
    try:
        function.visualize(grid, filename=tmpfile[:-4], subsampling=subsampling)
    except AttributeError:
        GridFunction(grid, function).visualize(grid, filename=tmpfile[:-4], subsampling=subsampling)
    return plot(tmpfile, color_attribute_name=function.name)
