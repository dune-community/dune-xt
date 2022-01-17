# ~~~
# This file is part of the dune-xt project:
#   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
# Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2017, 2019 - 2020)
#   René Fritze     (2018 - 2019)
#   Tim Keil        (2018, 2020)
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


    def visualize_function_on_dd_grid(function, dd_grid, subdomains=None):
        assert function.dim_domain == 2, f'Not implemented yet for {function.dim_domain}-dimensional grids!'
        assert function.dim_range == 1, f'Not implemented yet for {function.dim_domain}-dimensional functions!'
        subdomains = subdomains or list(range(dd_grid.num_subdomains))
        assert isinstance(subdomains, list), 'Please provide a list of subdomains'
        assert all(sd < dd_grid.num_subdomains for sd in subdomains)
        assert all(isinstance(sd, int) for sd in subdomains)
        tmpfile = tempfile.NamedTemporaryFile(mode='wb', delete=False, suffix='.vtu').name
        dd_grid.write_global_visualization(tmpfile[5:-4], function, subdomains)
        prestring = f's{dd_grid.num_subdomains:04d}-'
        if len(subdomains) == 1:
            prestring += f'p{subdomains[0]:04d}-'
            return plot(prestring + tmpfile[5:], color_attribute_name=function.name)
        try:
            plot(prestring + tmpfile[5:-4] + '.pvtu', color_attribute_name=function.name)
        except AttributeError:
            if 1 < len(subdomains) < dd_grid.num_subdomains:
                # TODO: fix this error
                print("For 1 < len(subdomains) < num_subdomains, the generated .pvtu file is known "
                       + "to be broken. The result for each subdomain can be viewed via paraview")
            else:
                print("Unexpected error")
