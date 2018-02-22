# This file is part of the dune-xt-la project:
#   https://github.com/dune-community/dune-xt-la
# Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Rene Milk (2017 - 2018)

from dune.xt import codegen


def matrices(cache):
  mat = ['CommonDenseMatrix', 'CommonSparseMatrixCsr', 'CommonSparseMatrixCsc',
    'CommonSparseOrDenseMatrixCsr', 'CommonSparseOrDenseMatrixCsc']

  if codegen.have_eigen(cache):
      mat.extend(['EigenRowMajorSparseMatrix', 'EigenDenseMatrix'])
  if codegen.have_istl(cache):
      mat.append('IstlRowMajorSparseMatrix')
  return mat


def vectors(cache):
  vecs = ['CommonDenseVector', 'CommonSparseVector']
  if codegen.have_eigen(cache):
      vecs.extend(['EigenDenseVector', 'EigenMappedDenseVector'])
  if codegen.have_istl(cache):
      vecs.append('IstlDenseVector')
  return vecs


def name_type_tuple(c, f):
    v = latype(c, f)
    return codegen.typeid_to_typedef_name(c+f), v


def latype(c, f):
    return 'Dune::XT::LA::{}<{}>'.format(c, f)


def commontype(c, f):
    return 'Dune::XT::Common::{}<{}>'.format(c, f)


def dunetype(c, f):
    return 'Dune::{}<{}>'.format(c, f)


def fieldtypes(cache):
    return ['double', 'std::complex<double>']


def vector_filter(vector, field):
    return not ('EigenMappedDenseVector' in vector  and 'complex' in field)
