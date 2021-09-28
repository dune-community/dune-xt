# ~~~
# This file is part of the dune-xt project:
#   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
# Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   René Fritze    (2019)
#   Tobias Leibner (2019 - 2020)
# ~~~

from itertools import product
from matrices import matrices, latype, fieldtypes
from dune.xt.codegen import typeid_to_typedef_name as safe_name

testtypes = [(safe_name('{}_{}'.format(m, f)), latype(m, f)) for m, f in product(matrices(cache), fieldtypes(cache))]
