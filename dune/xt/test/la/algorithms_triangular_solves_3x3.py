# ~~~
# This file is part of the dune-xt project:
#   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
# Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2016)
#   René Fritze     (2016 - 2019)
#   Tobias Leibner  (2015 - 2018, 2020)
# ~~~

from algorithms import la_types, common_types, dune_types, testtypes_creator

testtypes = testtypes_creator(la_types, common_types, dune_types, cache)
