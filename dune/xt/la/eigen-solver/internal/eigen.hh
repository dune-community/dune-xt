// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2017 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler  (2017)

#ifndef DUNE_XT_LA_EIGEN_SOLVER_INTERNAL_EIGEN_HH
#define DUNE_XT_LA_EIGEN_SOLVER_INTERNAL_EIGEN_HH

#if HAVE_EIGEN
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#endif

#include <dune/common/dynvector.hh>
#include <dune/common/typetraits.hh>

#include <dune/xt/la/container/eigen.hh>
#include <dune/xt/la/exceptions.hh>


namespace Dune {
namespace XT {
namespace LA {
namespace internal {


#if HAVE_EIGEN


template <class S>
std::vector<std::complex<S>>
compute_all_eigenvalues_using_eigen(const ::Eigen::Matrix<S, ::Eigen::Dynamic, ::Eigen::Dynamic>& matrix)
{
  ::Eigen::EigenSolver<::Eigen::Matrix<S, ::Eigen::Dynamic, ::Eigen::Dynamic>> eigen_solver(
      matrix, /*computeEigenvectors=*/false);
  if (eigen_solver.info() != ::Eigen::Success)
    DUNE_THROW(Exceptions::eigen_solver_failed, "The eigen backend reported '" << eigen_solver.info() << "'!");
  const auto evs = eigen_solver.eigenvalues(); // this should be an Eigen vector of std::complex<S>
  std::vector<std::complex<S>> ret(evs.size());
  for (size_t ii = 0; ii < size_t(evs.size()); ++ii)
    ret[ii] = evs[ii];
  return ret;
} // ... compute_all_eigenvalues_using_eigen(...)


template <class S>
std::vector<EigenDenseVector<std::complex<typename EigenDenseVector<S>::RealType>>>
compute_all_eigenvectors_using_eigen(const ::Eigen::Matrix<S, ::Eigen::Dynamic, ::Eigen::Dynamic>& matrix)
{
  ::Eigen::EigenSolver<::Eigen::Matrix<S, ::Eigen::Dynamic, ::Eigen::Dynamic>> eigen_solver(
      matrix, /*computeEigenvectors=*/true);
  if (eigen_solver.info() != ::Eigen::Success)
    DUNE_THROW(Exceptions::eigen_solver_failed, "The eigen backend reported '" << eigen_solver.info() << "'!");
  const auto eigenvectors = eigen_solver.eigenvectors(); // This should be an Eigen Matrix of std::complex<S> which
  //                                                        contains the eigenvectors in the columns.
  std::vector<EigenDenseVector<std::complex<typename EigenDenseVector<S>::RealType>>> ret;
  for (size_t ii = 0; ii < size_t(eigenvectors.cols()); ++ii)
    ret.emplace_back(eigenvectors.col(ii));
  return ret;
} // ... compute_all_eigenvectors_using_eigen(...)


#else // HAVE_EIGEN


template <class S>
std::vector<std::complex<S>>
compute_all_eigenvalues_using_eigen(const ::Eigen::Matrix<S, ::Eigen::Dynamic, ::Eigen::Dynamic>& /*matrix*/)
{
  static_assert(AlwaysFalse<S>::value, "You are missing Eigen!");
}

template <class S>
std::vector<EigenDenseVector<std::complex<typename EigenDenseVector<S>::RealType>>>
compute_all_eigenvectors_using_eigen(const ::Eigen::Matrix<S, ::Eigen::Dynamic, ::Eigen::Dynamic>& /*matrix*/)
{
  static_assert(AlwaysFalse<S>::value, "You are missing Eigen!");
}


#endif // HAVE_EIGEN


} // namespace internal
} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_EIGEN_SOLVER_INTERNAL_EIGEN_HH
