// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2013 - 2014, 2016)
//   Rene Milk       (2015)

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

template <class ScalarType, ChooseBackend backend = default_backend>
struct Container;

template <class ScalarType>
struct Container<ScalarType, ChooseBackend::common_dense>
{
  typedef CommonDenseVector<ScalarType> VectorType;
  typedef CommonDenseMatrix<ScalarType> MatrixType;
}; // struct Container< ..., common_dense >

template <class ScalarType>
struct Container<ScalarType, ChooseBackend::eigen_dense>
{
  typedef EigenDenseVector<ScalarType> VectorType;
  typedef EigenDenseMatrix<ScalarType> MatrixType;
}; // struct Container< ..., eigen_dense >

template <class ScalarType>
struct Container<ScalarType, ChooseBackend::eigen_sparse>
{
  typedef EigenDenseVector<ScalarType> VectorType;
  typedef EigenRowMajorSparseMatrix<ScalarType> MatrixType;
}; // struct Container< ..., eigen_sparse >

template <class ScalarType>
struct Container<ScalarType, ChooseBackend::istl_sparse>
{
  typedef IstlDenseVector<ScalarType> VectorType;
  typedef IstlRowMajorSparseMatrix<ScalarType> MatrixType;
}; // struct Container< ..., istl_sparse >

} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_CONTAINER_HH
