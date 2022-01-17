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
from dune.xt import guarded_import
from dune.xt.common.config import config

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


if config.HAVE_K3D:
    import os
    import tempfile
    from dune.xt.common.vtk.plot import plot
    from dune.xt.functions._functions_gridfunction import GridFunction

    def visualize_function(functions, grid, interactive=False, subsampling=False):
        if not isinstance(functions, (list, tuple)):
            functions = (functions,)
        # this is a (bad) implicit check that these are functions
        dim_domain = functions[0].dim_domain
        dim_range = functions[0].dim_range
        name = functions[0].name
        assert(all(f.dim_domain == dim_domain for f in functions))
        assert(all(f.dim_range == dim_range for f in functions))
        assert(all(f.name == name for f in functions))

        def visualize_single(func, filename):
            try:
                func.visualize(grid, filename=filename, subsampling=subsampling)
            except AttributeError:
                GridFunction(grid, func).visualize(grid, filename=filename, subsampling=subsampling)

        prefix = os.path.join(tempfile.gettempdir(), next(tempfile._get_candidate_names()))
        suffix = 'vt{}'.format('p' if dim_domain == 1 else 'u')
        if len(functions) == 0:
            return
        for ii, func in enumerate(functions):
            visualize_single(func, f'{prefix}_{ii}')
        if len(functions) == 1:
            filename = f'{prefix}_0.{suffix}'
        else:
            with open(f'{prefix}.pvd', 'w') as pvd_file:
                pvd_file.write('<?xml version=\'1.0\'?>\n')
                pvd_file.write('<VTKFile type=\'Collection\' version=\'0.1\'>\n')
                pvd_file.write('<Collection>\n')
                for ii, func in enumerate(functions):
                    pvd_file.write(
            f'<DataSet timestep=\'{ii}\' part=\'1\' name=\'{name}\' file=\'{prefix}_{ii}.{suffix}\'/>\n')
                pvd_file.write('</Collection>\n')
                pvd_file.write('</VTKFile>\n')
            filename = f'{prefix}.pvd'

        return plot(filename, color_attribute_name=name, interactive=True)

