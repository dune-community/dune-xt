// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2014, 2016 - 2017)
//   Rene Milk       (2015 - 2016, 2018)
//   Tobias Leibner  (2017)

#ifndef DUNE_XT_LA_CONTAINER_HH
#define DUNE_XT_LA_CONTAINER_HH

#include "container/interfaces.hh"
#include "container/common.hh"
#include "container/eigen.hh"
#include "container/istl.hh"

#include <dune/xt/common/logging.hh>

namespace Dune {
namespace XT {
namespace LA {


template <class ScalarType, Backends backend = default_backend>
struct Container;

template <class ScalarType>
struct Container<ScalarType, Backends::common_dense>
{
  typedef CommonDenseVector<ScalarType> VectorType;
  typedef CommonDenseMatrix<ScalarType> MatrixType;
}; // struct Container<..., common_dense>

template <class ScalarType>
struct Container<ScalarType, Backends::common_sparse>
{
  typedef CommonDenseVector<ScalarType> VectorType;
  typedef CommonSparseMatrix<ScalarType> MatrixType;
}; // struct Container<..., common_sparse>

template <class ScalarType>
struct Container<ScalarType, Backends::eigen_dense>
{
  typedef EigenDenseVector<ScalarType> VectorType;
  typedef EigenDenseMatrix<ScalarType> MatrixType;
}; // struct Container<..., eigen_dense>

template <class ScalarType>
struct Container<ScalarType, Backends::eigen_sparse>
{
  typedef EigenDenseVector<ScalarType> VectorType;
  typedef EigenRowMajorSparseMatrix<ScalarType> MatrixType;
}; // struct Container<..., eigen_sparse>

template <class ScalarType>
struct Container<ScalarType, Backends::istl_dense>
{
  typedef IstlDenseVector<ScalarType> VectorType;
  typedef IstlRowMajorSparseMatrix<ScalarType> MatrixType;
}; // struct Container<..., istl_dense>

template <class ScalarType>
struct Container<ScalarType, Backends::istl_sparse>
{
  typedef IstlDenseVector<ScalarType> VectorType;
  typedef IstlRowMajorSparseMatrix<ScalarType> MatrixType;
}; // struct Container<..., istl_sparse>


} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_CONTAINER_HH
