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

from dune.xt.grid import Dim, Cube, Simplex, make_cube_grid


init_args = (
        (Dim(1), [0], [1], [2]),
        (Dim(2), Cube(), [0, 0], [1, 1], [2, 2]),
        (Dim(3), Cube(), [0, 0, 0], [1, 1, 1], [2, 2, 2]),
        )


def test_init():
    for args in init_args:
        grid = make_cube_grid(*args)
