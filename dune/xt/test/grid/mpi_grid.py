# ~~~
# This file is part of the dune-xt project:
#   https://github.com/dune-community/dune-xt
# Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   René Fritze    (2018 - 2019)
#   Tobias Leibner (2020)
# ~~~

import dune.xt.grid.types as grid_types
from dune.xt.codegen import typeid_to_typedef_name as safe_name

# alberta needs manual flag adding in cmake, so we skip it here
all_grids = ((safe_name(g), g) for g in grid_types.all_types(cache, list(range(1, 4))) if 'Alberta' not in g)
