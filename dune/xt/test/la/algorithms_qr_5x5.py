# ~~~
# This file is part of the dune-xt project:
#   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
# Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   René Fritze    (2018 - 2019)
#   Tobias Leibner (2018 - 2020)
# ~~~

from algorithms import la_types, testtypes_creator

common_types = [f.split('_') for f in ['FieldMatrix_FieldVector_FieldVector_double,5,5']]

dune_types = [
    f.split('_')
    for f in ['FieldMatrix_FieldVector_FieldVector_double,5,5', 'DynamicMatrix_DynamicVector_DynamicVector_double']
]

testtypes = testtypes_creator(la_types, common_types, dune_types, cache)
