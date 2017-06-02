// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2017 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler  (2017)

#ifndef DUNE_XT_LA_EIGEN_SOLVER_EIGEN_HH
#define DUNE_XT_LA_EIGEN_SOLVER_EIGEN_HH

#if HAVE_EIGEN
#include <Eigen/Eigenvalues>
#endif

#include <dune/common/typetraits.hh>

#include <dune/xt/la/container/eigen/dense.hh>

#include "../eigen-solver.hh"

namespace Dune {
namespace XT {
namespace LA {

#if HAVE_EIGEN


template <class S>
class EigenSolver<EigenDenseMatrix<S>>
{
public:
  EigenSolver(const EigenDenseMatrix<S>& matrix)
    : matrix_(matrix)
  {
  }

  DynamicVector<S> all_eigenvalues() const
  {
    ::Eigen::EigenSolver<typename EigenDenseMatrix<S>::BackendType> eigen_solver(matrix_.backend());
    if (eigen_solver.info() != ::Eigen::Success)
      DUNE_THROW(Exceptions::eigen_solver_failed, "The eigen backend reported '" << eigen_solver.info() << "'!");
    const auto eigenvalues = eigen_solver.eigenvalues(); // <- this should be an Eigen vector of std::complex
    DynamicVector<S> ret(eigenvalues.size());
    for (size_t ii = 0; ii < eigenvalues.size(); ++ii) {
      if (std::abs(eigenvalues[ii].imag()) > 1e-15)
        DUNE_THROW(
            Exceptions::eigen_solver_failed,
            "Not all eigenvalues are purely real!\n   eigenvalues[" << ii << "].imag(): " << eigenvalues[ii].imag());
      if (eigenvalues[ii].real() < 1e-15)
        DUNE_THROW(
            Exceptions::eigen_solver_failed,
            "Not all eigenvalues are positive!\n   eigenvalues[" << ii << "].real(): " << eigenvalues[ii].real());
      ret[ii] = eigenvalues[ii].real();
    }
    return ret;
  }

  S min_eigenvalue() const
  {
    S min_ev = std::numeric_limits<S>::max();
    for (const auto& ev : all_eigenvalues())
      min_ev = std::min(ev, min_ev);
    return min_ev;
  }

  S max_eigenvalue() const
  {
    S max_ev = std::numeric_limits<S>::min();
    for (const auto& ev : all_eigenvalues())
      max_ev = std::max(ev, max_ev);
    return max_ev;
  }

private:
  const EigenDenseMatrix<S>& matrix_;
}; // class EigenSolver<EigenDenseMatrix<S>>


#else // HAVE_EIGEN


template <class S>
class EigenSolver<EigenDenseMatrix<S>>
{
  static_assert(AlwaysFalse<S>::value, "You are missing eigen!")
};


#endif // HAVE_EIGEN

} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_EIGEN_SOLVER_EIGEN_HH
