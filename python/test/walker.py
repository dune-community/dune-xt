# ~~~
# This file is part of the dune-xt-grid project:
#   https://github.com/dune-community/dune-xt-grid
# Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   René Fritze    (2018)
#   Tobias Leibner (2018)
# ~~~

import functools
import itertools

import pytest
import dune.xt.common as xtc
import dune.xt.grid as xtg
from dune.xt.grid import provider, types, walker


def _grid_provider_factory(grid_type, mpi, maker_str='make_cube_grid__{}'):
    if not mpi:
        fn = maker_str.format(grid_type)
        maker = getattr(xtg.provider, fn)
        return maker()

    try:
        from mpi4py import MPI
    except ImportError:
        pytest.skip('optional mpi4py is missing')
        return
    opts = provider.default_options_cube_grid(grid_type)
    fn = maker_str.format(grid_type)
    maker = getattr(xtg.provider, fn)
    return maker(opts, MPI.COMM_WORLD)


_dd_subdomain_grid_provider_factory = functools.partial(
    _grid_provider_factory, maker_str='make_cube_dd_subdomain_grid__{}')
'''
@pytest.fixture(params=xtg.types.available_types)
def mpi_grid_provider(request):
    return _grid_provider_factory(request.param, mpi=True)


@pytest.fixture(params=xtg.types.available_types)
def grid_provider(request):
    return _grid_provider_factory(request.param, mpi=False)


@pytest.fixture(params=itertools.product(xtg.types.available_types, (True, False)))
def combined_grid_provider(request):
    return _grid_provider_factory(*request.param)


@pytest.fixture(params=itertools.product(xtg.types.available_types, (True, False)))
def combined_dd_subdomain_grid_provider(request):
    return _grid_provider_factory(*request.param)
'''


def test_available():
    assert len(xtg.types.available_types) > 0


'''
def test_grid_provider(combined_grid_provider):
    grid_provider = combined_grid_provider
    assert grid_provider.max_level() >= 0
    num_el = grid_provider.num_elements
    assert num_el > 1
    grid_provider.global_refine(1)
    try:
        grid_provider.visualize()
    except xtc.DuneError as e:
        if 'NotImplemented' not in str(e):
            raise e
    assert grid_provider.num_elements > num_el
    assert grid_provider.num_subdomains == 1
    assert grid_provider.max_entity_diameter() > 0
    assert grid_provider.dim > 0


def test_dd_subdomain_grid_provider(combined_dd_subdomain_grid_provider):
    grid_provider = combined_dd_subdomain_grid_provider
    assert grid_provider.max_level() >= 0
    num_el = grid_provider.num_elements
    assert num_el > 1
    grid_provider.global_refine(1)
    try:
        grid_provider.visualize()
    except xtc.DuneError as e:
        if 'NotImplemented' not in str(e):
            raise e
    assert grid_provider.num_elements > num_el
    assert grid_provider.num_subdomains == 1
    assert grid_provider.max_entity_diameter() > 0
    assert grid_provider.dim > 0


def test_walker(grid_provider):
    walker = xtg.walker.make_walker(grid_provider)
    walker.walk()
    walker.clear()
'''


def test_count():
    pass
