// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   René Fritze     (2018 - 2019)
//   Tobias Leibner  (2018, 2020)

#ifndef DUNE_XT_LA_MATRIX_INVERTER_INTERNAL_EIGEN_HH
#define DUNE_XT_LA_MATRIX_INVERTER_INTERNAL_EIGEN_HH

#include <limits>

#if HAVE_EIGEN
#  include <dune/xt/common/disable_warnings.hh>
#  include <Eigen/Core>
#  include <Eigen/SVD>
// needed for the inverse method of Eigen::Matrix, see http://eigen.tuxfamily.org/bz/show_bug.cgi?id=1089
#  include <Eigen/LU>
#  include <dune/xt/common/reenable_warnings.hh>
#endif

#include <dune/common/typetraits.hh>

namespace Dune::XT::LA::internal {

#if HAVE_EIGEN


/**
 * \brief Computes the Moore-Penrose inverse.
 * \sa    https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_inverse
 *
 *        The implementation is inspired by \sa http://eigen.tuxfamily.org/bz/show_bug.cgi?id=257#c14
 */
template <class S>
::Eigen::Matrix<S, ::Eigen::Dynamic, ::Eigen::Dynamic>
compute_moore_penrose_inverse_using_eigen(const ::Eigen::Matrix<S, ::Eigen::Dynamic, ::Eigen::Dynamic>& matrix,
                                          const double& epsilon = std::numeric_limits<double>::epsilon())
{
  ::Eigen::JacobiSVD<::Eigen::Matrix<S, ::Eigen::Dynamic, ::Eigen::Dynamic>> svd(
      matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
  const double tolerance = epsilon * std::max(matrix.cols(), matrix.rows()) * svd.singularValues().array().abs()(0);
  return svd.matrixV()
         * (svd.singularValues().array().abs() > tolerance)
               .select(svd.singularValues().array().inverse(), 0)
               .matrix()
               .asDiagonal()
         * svd.matrixU().adjoint();
}


#else // HAVE_EIGEN


template <class M>
M compute_moore_penrose_inverse_using_eigen(const M& matrix,
                                            const double& /*epsilon*/ = std::numeric_limits<double>::epsilon())
{
  static_assert(AlwaysFalse<M>::value, "You are missing eigen!");
  return matrix;
}


#endif // HAVE_EIGEN

} // namespace Dune::XT::LA::internal

#endif // DUNE_XT_LA_MATRIX_INVERTER_INTERNAL_EIGEN_HH
