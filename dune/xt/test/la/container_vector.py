# ~~~
# This file is part of the dune-xt project:
#   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
# Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2016)
#   René Fritze     (2017 - 2019)
#   Tobias Leibner  (2015, 2019 - 2020)
# ~~~

from itertools import product
from matrices import latype, vectors, fieldtypes, vector_filter
from dune.xt.codegen import typeid_to_typedef_name as safe_name

testtypes = [(safe_name('{}_{}'.format(mv, f)), latype(mv, f))
             for mv, f in product(vectors(cache), fieldtypes(cache))
             if vector_filter(mv, f)]
