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
    from dune.xt.grid._boundaryinfo import *
except ImportError as e:
    import os
    import logging
    if os.environ.get('DXT_PYTHON_DEBUG', False):
        raise e
    logging.error('dune-xt-grid bindings not available')


def _meta_boundary(name, grid_provider, config):
    for factory in [globals()[s] for s in globals().keys() if s.startswith(name)]:
        try:
            return factory(grid_provider, config)
        except:
            continue
    raise TypeError('no matching {} for boundaryinfo {}'.format(name, grid_provider.__class__))


def make_boundary_info_on_dd_subdomain_boundary_layer(grid_provider, config):
    return _meta_boundary('make_boundary_info_on_dd_subdomain_boundary_layer', grid_provider, config)


def make_boundary_info_on_dd_subdomain_layer(grid_provider, config):
    return _meta_boundary('make_boundary_info_on_dd_subdomain_layer', grid_provider, config)


def make_boundary_info_on_leaf_layer(grid_provider, config):
    return _meta_boundary('make_boundary_info_on_leaf_layer', grid_provider, config)
