# ~~~
# This file is part of the dune-xt project:
#   https://github.com/dune-community/dune-xt
# Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   René Fritze    (2018 - 2019)
#   Tobias Leibner (2019 - 2020)
# ~~~

try:
    from dune.xt.grid._provider import *
except ImportError as e:
    import os
    import logging
    if os.environ.get('DXT_PYTHON_DEBUG', False):
        raise e
    logging.error('dune-xt-grid bindings not available')


def default_options_cube_grid(type_str):
    factory = 'GridProviderFactory__{}'.format(type_str)
    try:
        fac = globals()[factory]
    except KeyError:
        raise TypeError('no GridProviderFactory available for Type {}'.format(type_str))
    return fac.default_config('xt.grid.gridprovider.cube')
