# ~~~
# This file is part of the dune-xt project:
#   https://github.com/dune-community/dune-xt
# Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2016 - 2017)
#   René Fritze     (2018 - 2019)
#   Tim Keil        (2018)
#   Tobias Leibner  (2018 - 2020)
# ~~~

import pytest

from dune.xt.grid import Dim, Cube, make_cube_grid
from dune.xt.test._test_grid_intersection import call_on_each_intersection


def test_print():
    grid = make_cube_grid(Dim(2), Cube(), [0, 0], [1, 1], [2, 2])
    call_on_each_intersection(grid, lambda intersection: print(intersection))


def test_boundary():
    grid = make_cube_grid(Dim(2), Cube(), [0, 0], [1, 1], [2, 2])
    call_on_each_intersection(grid, lambda intersection: print(intersection.boundary))


def test_boundary_segment_index():
    grid = make_cube_grid(Dim(2), Cube(), [0, 0], [1, 1], [2, 2])
    call_on_each_intersection(
        grid, lambda intersection: intersection.boundary and print(intersection.boundary_segment_index))


def test_neighbor():
    grid = make_cube_grid(Dim(2), Cube(), [0, 0], [1, 1], [2, 2])
    call_on_each_intersection(grid, lambda intersection: print(intersection.neighbor))


def test_index_in_inside():
    grid = make_cube_grid(Dim(2), Cube(), [0, 0], [1, 1], [2, 2])
    call_on_each_intersection(grid, lambda intersection: print(intersection.index_in_inside))


def test_index_in_outside():
    grid = make_cube_grid(Dim(2), Cube(), [0, 0], [1, 1], [2, 2])
    call_on_each_intersection(grid, lambda intersection: print(intersection.index_in_outside))


def test_center_unit_outer_normal():
    grid = make_cube_grid(Dim(2), Cube(), [0, 0], [1, 1], [2, 2])
    call_on_each_intersection(grid, lambda intersection: print(intersection.center_unit_outer_normal))
