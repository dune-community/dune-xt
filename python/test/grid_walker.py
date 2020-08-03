# ~~~
# This file is part of the dune-xt project:
#   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
# Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2020)
#   René Fritze     (2018 - 2019)
#   Tobias Leibner  (2018 - 2020)
# ~~~

from dune.xt.test.base import runmodule
from dune.xt.grid import Dim, Cube, Simplex, make_cube_grid, Walker

grid_elements_per_dim = 10

init_args = (
    (Dim(1), [0], [1], [grid_elements_per_dim]),
    (Dim(2), Cube(), [0, 0], [1, 1], [grid_elements_per_dim]*2),
    (Dim(2), Simplex(), [0, 0], [1, 1], [grid_elements_per_dim]*2),
    (Dim(3), Cube(), [0, 0, 0], [1, 1, 1], [grid_elements_per_dim]*3),
    (Dim(3), Simplex(), [0, 0, 0], [1, 1, 1], [grid_elements_per_dim]*3),
)


def test_init():
    for args in init_args:
        grid = make_cube_grid(*args)
        walker = Walker(grid)


def test_walk():
    for args in init_args:
        grid = make_cube_grid(*args)
        walker = Walker(grid)
        walker.walk(True)


if __name__ == "__main__":
    runmodule(filename=__file__)
