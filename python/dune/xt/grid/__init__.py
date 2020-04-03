# ~~~
# This file is part of the dune-xt project:
#   https://github.com/dune-community/dune-xt
# Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2017, 2019)
#   René Fritze     (2018 - 2019)
#   Tobias Leibner  (2019 - 2020)
# ~~~

from dune.xt import guarded_import

for mod_name in (
        '_boundaryinfo',
        '_types',
        '_walker',
        '_provider',
):
    guarded_import(globals(), 'dune.xt.grid', mod_name)


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
