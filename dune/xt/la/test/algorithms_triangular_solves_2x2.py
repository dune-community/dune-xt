# ~~~
# This file is part of the dune-xt-la project:
#   https://github.com/dune-community/dune-xt-la
# Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Tobias Leibner (2018)
# ~~~

from algorithms import la_types, testtypes_creator

common_types = [f.split('_') for f in ['FieldMatrix_FieldVector_FieldVector_double,2,2']]

dune_types = [f.split('_') for f in ['FieldMatrix_FieldVector_FieldVector_double,2,2',
                                     'DynamicMatrix_DynamicVector_DynamicVector_double']]

testtypes = testtypes_creator(la_types, common_types, dune_types, cache)
