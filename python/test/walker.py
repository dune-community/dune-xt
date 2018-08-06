# ~~~
# This file is part of the dune-xt-grid project:
#   https://github.com/dune-community/dune-xt-grid
# Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Rene Milk (2018)
# ~~~

import pytest
import dune.xt.common as xtc
import dune.xt.grid as xtg

@pytest.fixture(params=xtg.available_types)
def mpi_grid_provider(request):
    try:
        from mpi4py import MPI
    except ImportError:
        pytest.skip('optional mpi4py is missing')
        return
    fn = 'make_cube_grid__{}'.format(request.param)
    maker = getattr(xtg, fn)
    return maker()

@pytest.fixture(params=xtg.available_types)
def grid_provider(request):
    fn = 'make_cube_grid__{}'.format(request.param)
    maker = getattr(xtg, fn)
    return maker()


def test_available():
    assert len(xtg.available_types) > 0


def test_grid_provider(grid_provider):
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
    walker = xtg.make_walker(grid_provider)
    walker.walk()
    walker.clear()



def test_mpi_cubegrid(mpi_grid_provider):
    pass

def test_count():
    pass

