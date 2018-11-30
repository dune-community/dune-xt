# ~~~
# This file is part of the dune-xt-la project:
#   https://github.com/dune-community/dune-xt-la
# Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   René Fritze (2017 - 2018)
# ~~~

from itertools import product

from dune.xt import codegen

import matrices


conts = matrices.matrices(cache) + matrices.vectors(cache)
container = [matrices.name_type_tuple(c,f)
             for c,f in product(conts, matrices.fieldtypes(cache))
             if matrices.vector_filter(c,f)]
