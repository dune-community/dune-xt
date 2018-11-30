# ~~~
# This file is part of the dune-xt-la project:
#   https://github.com/dune-community/dune-xt-la
# Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2016)
#   René Fritze     (2017 - 2018)
#   Tobias Leibner  (2015)
# ~~~

from itertools import product
from matrices import matrices, latype, vectors, fieldtypes, vector_filter
from dune.xt.codegen import typeid_to_typedef_name as safe_name

testtypes = [(safe_name('{}_{}_{}'.format(*mv,f)), latype(mv[0],f), latype(mv[1],f))
             for mv,f in product(zip(matrices(cache), vectors(cache)), fieldtypes(cache))
             if vector_filter(mv[1], f)]