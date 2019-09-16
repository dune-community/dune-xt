# ~~~
# This file is part of the dune-xt project:
#   https://github.com/dune-community/dune-xt
# Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   René Fritze    (2018 - 2019)
#   Tim Keil       (2018)
#   Tobias Leibner (2018 - 2019)
# ~~~

import pytest
import dune.xt.grid as xtg
from dune.xt.grid import types
import dune.xt.functions as xtf


@pytest.fixture(params=xtg.types.available_types)
def function_provider(request):
    grid_type = request.param
    dim = 1
    if '2d' in grid_type:
        dim = 2
    if '3d' in grid_type:
        dim = 3
    fn = "ConstantFunction__{}d_to_{}x{}".format(dim, 1, 1)
    maker = getattr(xtf, fn)
    return maker([2], 'test_function')


def test_available():
    assert len(xtg.types.available_types) > 0


def test_grid_provider(function_provider):
    assert function_provider.name == 'test_function'


def test_count():
    pass
