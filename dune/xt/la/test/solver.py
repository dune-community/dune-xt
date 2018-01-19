# This file is part of the dune-xt-la project:
#   https://github.com/dune-community/dune-xt-la
# Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2016)
#   Rene Milk       (2017 - 2018)
#   Tobias Leibner  (2015 - 2017)

from matrices import latype
from dune.xt.codegen import have_eigen, have_istl, typeid_to_typedef_name as safe_name


types = [f.split('_') for f in ['CommonDenseMatrix_CommonDenseVector_CommonDenseVector_complex',
                        'CommonDenseMatrix_CommonDenseVector_CommonDenseVector_double',
                        'EigenDenseMatrix_EigenDenseVector_EigenDenseVector_complex',
                        'EigenDenseMatrix_EigenDenseVector_EigenDenseVector_double',
                        'EigenDenseMatrix_EigenDenseVector_EigenMappedDenseVector_double',
                        'EigenDenseMatrix_EigenMappedDenseVector_EigenDenseVector_double',
                        'EigenDenseMatrix_EigenMappedDenseVector_EigenMappedDenseVector_double',
                        'EigenRowMajorSparseMatrix_EigenDenseVector_EigenDenseVector_complex',
                        'EigenRowMajorSparseMatrix_EigenDenseVector_EigenDenseVector_double',
                        'IstlRowMajorSparseMatrix_IstlDenseVector_IstlDenseVector_double']]


def test_tuple(args):
    o, r, s, f = args
    if f == 'complex':
        f = 'std::complex<double>'
    return (safe_name('{}_{}_{}_{}'.format(o, r, s, f)), latype(o,f), latype(r,f), latype(s,f))


def type_ok(t):
    if sum(['Eigen' in x for x in t]):
        return have_eigen(cache)
    if sum(['Istl' in x for x in t]):
        return have_istl(cache)
    return True


testtypes = [test_tuple(item) for item in types if type_ok(item)]
