// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_LA_CONTAINER_HH
#define DUNE_STUFF_LA_CONTAINER_HH

#include <dune/common/typetraits.hh>

#include "container/interfaces.hh"
#include "container/common.hh"
#include "container/eigen.hh"
#include "container/istl.hh"

namespace Dune {
namespace Stuff {
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
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_LA_CONTAINER_HH
