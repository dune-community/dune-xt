# This file is part of the dune-xt-la project:
#   https://github.com/dune-community/dune-xt-la
# Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2016)
#   Rene Milk       (2017 - 2018)
#   Tobias Leibner  (2015 - 2017)

from algorithms import la_types, common_types, dune_types, testtypes_creator

testtypes = testtypes_creator(la_types, common_types, dune_types, cache)
