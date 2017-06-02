// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2017 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_XT_LA_EIGEN_SOLVER_HH
#define DUNE_XT_LA_EIGEN_SOLVER_HH

#include <dune/common/dynvector.hh>
#include <dune/common/typetraits.hh>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/la/type_traits.hh>

namespace Dune {
namespace XT {
namespace LA {
namespace Exceptions {


class eigen_solver_failed : public Dune::Exception
{
};


} // namespace Exceptions


template <class MatrixImp>
class EigenSolver
{
  static_assert(
      AlwaysFalse<MatrixImp>::value,
      "This is the unspecialized version of LA::EigenSolver< ... >, please add a specialization for this matrix!");

public:
  typedef MatrixImp MatrixType;
  typedef typename MatrixType::ScalarType ScalarType;

  EigenSolver(const MatrixType& /*matrix*/)
  {
    DUNE_THROW(NotImplemented,
               "This is the unspecialized version of LA::EigenSolver< ... >, please add a specialization for '"
                   << Common::Typename<MatrixType>::value()
                   << "'!");
  }

  DynamicVector<ScalarType> all_eigenvalues() const
  {
    DUNE_THROW(NotImplemented,
               "This is the unspecialized version of LA::EigenSolver< ... >, please add a specialization for '"
                   << Common::Typename<MatrixType>::value()
                   << "'!");
    return DenseVector<ScalarType>();
  }

  ScalarType min_eigenvalue() const
  {
    DUNE_THROW(NotImplemented,
               "This is the unspecialized version of LA::EigenSolver< ... >, please add a specialization for '"
                   << Common::Typename<MatrixType>::value()
                   << "'!");
    return ScalarType();
  }

  ScalarType max_eigenvalue() const
  {
    DUNE_THROW(NotImplemented,
               "This is the unspecialized version of LA::EigenSolver< ... >, please add a specialization for '"
                   << Common::Typename<MatrixType>::value()
                   << "'!");
    return ScalarType();
  }
}; // class EigenSolver


template <class M>
typename std::enable_if<XT::LA::is_matrix<M>::value, EigenSolver<M>>::type make_eigen_solver(const M& matrix)
{
  return EigenSolver<M>(matrix);
}


} // namespace LA
} // namespace XT
} // namespace Dune

#include "eigen-solver/eigen.hh"

#endif // DUNE_XT_LA_EIGEN_SOLVER_HH
