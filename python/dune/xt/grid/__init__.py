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
#   Tobias Leibner  (2019 - 2020)
# ~~~

from numbers import Number
import os
from tempfile import NamedTemporaryFile

from dune.xt import guarded_import
from dune.xt.common.vtk.plot import plot


for mod_name in (
        '_grid_boundaryinfo_alldirichlet',
        '_grid_boundaryinfo_allneumann',
        '_grid_boundaryinfo_allreflecting',
        '_grid_boundaryinfo_interfaces',
        '_grid_boundaryinfo_normalbased',
        '_grid_boundaryinfo_types',
        '_grid_filters_base',
        '_grid_filters_element',
        '_grid_filters_intersection',
        '_grid_functors_boundary_detector',
     # '_grid_functors_bounding_box',
        '_grid_functors_interfaces',
     # '_grid_functors_refinement',
        '_grid_gridprovider_cube',
        '_grid_gridprovider_gmsh',
        '_grid_gridprovider_provider',
        '_grid_intersection',
        '_grid_traits',
        '_grid_walker',
):
    guarded_import(globals(), 'dune.xt.grid', mod_name)


def Dim(d):
    assert isinstance(d, Number)
    if f'Dimension{d}' not in globals():
        raise RuntimeError(f'Dimension {d} not available, extend <python/dune/xt/grid/traits.cc>!')
    return globals()[f'Dimension{d}']()


def visualize_grid(grid):
    tmpfile = NamedTemporaryFile(mode='wb', delete=False, suffix='.vtu').name
    grid.visualize(tmpfile[:-4])
    return plot(tmpfile, color_attribute_name='Element index') # see visualize in python/dune/xt/grid/gridprovider.hh

