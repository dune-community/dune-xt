// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_LA_CONTAINER_HH
#define DUNE_STUFF_LA_CONTAINER_HH

#include "container/interfaces.hh"
#include "container/common.hh"
#include "container/eigen.hh"
#include "container/istl.hh"

#include <dune/stuff/common/logging.hh>

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

// template< class Traits, size_t domainDim, size_t rangeDim, size_t rangeDimCols > class SpaceInterface
template <class Space>
typename Space::RangeFieldType
communicated_dot(const Dune::Stuff::LA::IstlDenseVector<typename Space::RangeFieldType>& vector,
                 const Dune::Stuff::LA::IstlDenseVector<typename Space::RangeFieldType>& source, const Space& space)
{
  typename Space::RangeFieldType result = typename Space::RangeFieldType(0);
  space.communicator().dot(vector.backend(), source.backend(), result);
  return result;
}

template <template <class> class VectorImp, class Space>
typename Space::RangeFieldType communicated_dot(const VectorImp<typename Space::RangeFieldType>& vector,
                                                const VectorImp<typename Space::RangeFieldType>& source,
                                                const Space& /*space*/)
{
  auto result = vector.dot(source);
  DSC_LOG_DEBUG_0 << "communicated_dot does not account for overlapping dofs with non-ISTL vector types atm\n";
  return MPIHelper::getCollectiveCommunication().sum(result);
}

template <class Space>
void ensure_parallel_dof_consistency(Dune::Stuff::LA::IstlDenseVector<typename Space::RangeFieldType>& vector,
                                     const Space& space)
{
  space.communicator().copyOwnerToAll(vector.backend(), vector.backend());
}

template <template <class> class VectorImp, class Space>
void ensure_parallel_dof_consistency(VectorImp<typename Space::RangeFieldType>& /*vector*/, const Space& /*space*/)
{
  DSC_LOG_DEBUG_0 << "parallel dof consistency can only be ensured for ISTL vectors atm\n";
}

} // namespace LA
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_LA_CONTAINER_HH
