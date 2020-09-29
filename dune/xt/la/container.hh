// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2014, 2016 - 2017, 2019)
//   René Fritze     (2015 - 2016, 2018 - 2019)
//   Tobias Leibner  (2017 - 2018, 2020)

#ifndef DUNE_XT_LA_CONTAINER_HH
#define DUNE_XT_LA_CONTAINER_HH

#include <boost/tuple/tuple.hpp>

#include "container/interfaces.hh"
#include "container/common.hh"
#include "container/eigen.hh"
#include "container/istl.hh"

namespace Dune {
namespace XT {
namespace LA {


template <class ScalarType, Backends backend = default_backend>
struct Container;

template <class ScalarType>
struct Container<ScalarType, Backends::common_dense>
{
  using VectorType = CommonDenseVector<ScalarType>;
  using MatrixType = CommonDenseMatrix<ScalarType>;
}; // struct Container<..., common_dense>

template <class ScalarType>
struct Container<ScalarType, Backends::common_sparse>
{
  using VectorType = CommonDenseVector<ScalarType>;
  using MatrixType = CommonSparseMatrix<ScalarType>;
}; // struct Container<..., common_sparse>

template <class ScalarType>
struct Container<ScalarType, Backends::eigen_dense>
{
  using VectorType = EigenDenseVector<ScalarType>;
  using MatrixType = EigenDenseMatrix<ScalarType>;
}; // struct Container<..., eigen_dense>

template <class ScalarType>
struct Container<ScalarType, Backends::eigen_sparse>
{
  using VectorType = EigenDenseVector<ScalarType>;
  using MatrixType = EigenRowMajorSparseMatrix<ScalarType>;
}; // struct Container<..., eigen_sparse>

template <class ScalarType>
struct Container<ScalarType, Backends::istl_dense>
{
  using VectorType = IstlDenseVector<ScalarType>;
  using MatrixType = IstlRowMajorSparseMatrix<ScalarType>;
}; // struct Container<..., istl_dense>

template <class ScalarType>
struct Container<ScalarType, Backends::istl_sparse>
{
  using VectorType = IstlDenseVector<ScalarType>;
  using MatrixType = IstlRowMajorSparseMatrix<ScalarType>;
}; // struct Container<..., istl_sparse>


template <class S>
using AvailableVectorTypes = std::tuple<CommonDenseVector<S>,
                                        IstlDenseVector<S>
#if HAVE_EIGEN
                                        ,
                                        EigenDenseVector<S>
#endif
                                        >;


template <class S>
using AvailableDenseMatrixTypes = std::tuple<CommonDenseMatrix<S>
#if HAVE_EIGEN
                                             ,
                                             EigenDenseMatrix<S>
#endif
                                             >;


template <class S>
using AvailableSparseMatrixTypes = std::tuple<CommonSparseMatrix<S>,
                                              IstlRowMajorSparseMatrix<S>
#if HAVE_EIGEN
                                              ,
                                              EigenRowMajorSparseMatrix<S>
#endif
                                              >;


} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_CONTAINER_HH
