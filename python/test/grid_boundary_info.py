# ~~~
# This file is part of the dune-xt project:
#   https://github.com/dune-community/dune-xt
# Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   René Fritze    (2018 - 2019)
#   Tobias Leibner (2018 - 2019)
# ~~~

import pytest


def test_types():
    from dune.xt.grid import (
            NoBoundary,
            UnknownBoundary,
            DirichletBoundary,
            NeumannBoundary,
            RobinBoundary,
            ReflectingBoundary,
            AbsorbingBoundary,
            InflowBoundary,
            OutflowBoundary,
            InflowOutflowBoundary,
            ImpermeableBoundary,
            )
    NoBoundary()
    UnknownBoundary()
    DirichletBoundary()
    NeumannBoundary()
    RobinBoundary()
    ReflectingBoundary()
    AbsorbingBoundary()
    InflowBoundary()
    OutflowBoundary()
    InflowOutflowBoundary()
    ImpermeableBoundary()


def test_initless_boundary_infos():
    from dune.xt.grid import (
            AllDirichletBoundaryInfo,
            AllNeumannBoundaryInfo,
            AllReflectingBoundaryInfo,
            )
    from dune.xt.grid import make_cube_grid
    from grid_provider_cube import init_args as grid_init_args
    for args in grid_init_args:
        grid = make_cube_grid(*args)
        AllDirichletBoundaryInfo(grid)
        AllNeumannBoundaryInfo(grid)
        AllReflectingBoundaryInfo(grid)


def test_normalbased_boundary_inf():
    from dune.xt.grid import NormalBaseBoundaryInfo
    from dune.xt.grid import make_cube_grid, NoBoundary
    from grid_provider_cube import init_args as grid_init_args
    for args in grid_init_args:
        grid = make_cube_grid(*args)
        NormalBaseBoundaryInfo(
                grid_provider=grid,
                default_boundary_type=NoBoundary(),
                tolerance=1e-10,
                logging_prefix='')
        NormalBaseBoundaryInfo(
                grid_provider=grid,
                default_boundary_type=NoBoundary(),
                tolerance=1e-10)
        NormalBaseBoundaryInfo(
                grid_provider=grid,
                default_boundary_type=NoBoundary())
