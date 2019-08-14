# ~~~
# This file is part of the dune-xt-functions project:
#   https://github.com/dune-community/dune-xt-functions
# Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   René Fritze (2018 - 2019)
# ~~~

import pytest
from dune.xt.common.test import load_all_submodule
from dune.xt.grid import types, provider
import dune.xt.grid as xtg
import dune.xt.functions as xtf


def test_load_all():
    import dune.xt.functions as xtc
    load_all_submodule(xtc)


# @pytest.fixture(params=xtg.types.available_types)
# def _grid_provider_factory(request):
#     maker_str='make_cube_grid__{}'
#     fn = maker_str.format(request.param)
#     maker = getattr(xtg.provider, fn)
#     return maker()
#
#
# def test_create_functions(_grid_provider_factory):
#     grid = _grid_provider_factory
#     func = xtf.make_expression_function_1x1('x', '1+(cos(0.5*pi*x[0])*cos(0.5*pi*x[1]))', order=2, name='lambda_0')
#     func = xtf.make_constant_function_1x1(1)
#     func = xtf.make_constant_function_2x1(1)
