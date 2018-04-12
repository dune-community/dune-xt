# ~~~
# This file is part of the dune-xt-la project:
#   https://github.com/dune-community/dune-xt-la
# Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Tobias Leibner (2018)
# ~~~

from dune.xt.codegen import typeid_to_typedef_name as safe_name, have_eigen

matrix = ['EigenDenseMatrix<double>', 'FieldMatrix<double, 6, 6>']

def _ok(ft):
    if 'Eigen' in ft[0]:
        return have_eigen(cache)
    return True

testtypes = [(safe_name('_'.join(ft)), ft) for ft in matrix if _ok(ft)]
