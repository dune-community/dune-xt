# ~~~
# This file is part of the dune-xt-grid project:
#   https://github.com/dune-community/dune-xt-grid
# Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   René Fritze (2018 - 2019)
# ~~~

try:
    from dune.xt._walker import *
except ImportError as e:
    import os
    import logging
    if os.environ.get('DXT_PYTHON_DEBUG', False):
        raise e
    logging.error('dune-xt-grid bindings not available')


def make_walker(gridprovider, level=0):
    for factory in [globals()[s] for s in globals().keys() if s.startswith('make_walker_on_')]:
        try:
            return factory(gridprovider, level)
        except:
            continue
    raise TypeError('no matching walker for gridview {}'.format(gridprovider.__class__))


def make_apply_on_dirichlet_intersections(boundaryinfo, grid, layer='leaf_view', *args, **kwargs):
    factory = globals()['make_apply_on_dirichlet_intersections_{}_{}'.format(layer, grid.grid_type)]
    return factory(boundaryinfo, *args, **kwargs)
