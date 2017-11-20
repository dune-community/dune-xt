// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2017 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_XT_LA_EIGEN_SOLVER_EIGEN_HH
#define DUNE_XT_LA_EIGEN_SOLVER_EIGEN_HH

#include <algorithm>
#include <complex>
#include <functional>

#include <dune/xt/la/container.hh>
#include <dune/xt/la/solver.hh>
#include <dune/xt/la/eigen-solver.hh>

#include "internal/base.hh"
#include "internal/eigen.hh"
#include "internal/lapacke.hh"

namespace Dune {
namespace XT {
namespace LA {

#if HAVE_EIGEN


template <class S>
class EigenSolverOptions<EigenDenseMatrix<S>>
{
public:
  static std::vector<std::string> types()
  {
    return
    {
      "eigen"
#if HAVE_LAPACKE
          ,
          "lapack"
#endif
    };
  }

  static Common::Configuration options(const std::string type = "")
  {
    const std::string actual_type = type.empty() ? types()[0] : type;
    internal::ensure_eigen_solver_type(actual_type, types());
    Common::Configuration opts = internal::default_eigen_solver_options();
    opts["type"] = actual_type;
    return opts;
  }
}; // class EigenSolverOptions<EigenDenseMatrix<S>>


template <class S>
class EigenSolver<EigenDenseMatrix<S>> : public internal::EigenSolverBase<EigenDenseMatrix<S>,
                                                                          S,
                                                                          EigenDenseMatrix<XT::Common::real_t<S>>,
                                                                          EigenDenseMatrix<XT::Common::complex_t<S>>>
{
  using BaseType = internal::EigenSolverBase<EigenDenseMatrix<S>,
                                             S,
                                             EigenDenseMatrix<XT::Common::real_t<S>>,
                                             EigenDenseMatrix<XT::Common::complex_t<S>>>;

public:
  using typename BaseType::RealType;

  template <class... Args>
  explicit EigenSolver(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }

protected:
  void compute() const override final
  {
    const auto type = options_.template get<std::string>("type");
    const size_t N = matrix_.rows();
    if (type == "eigen") {
      if (options_.template get<bool>("compute_eigenvalues") && options_.template get<bool>("compute_eigenvectors")) {
        eigenvalues_ = std::make_unique<std::vector<XT::Common::complex_t<RealType>>>(N);
        eigenvectors_ = std::make_unique<EigenDenseMatrix<XT::Common::complex_t<S>>>(N, N);
        internal::compute_all_eigenvalues_and_vectors_using_eigen(
            matrix_.backend(), *eigenvalues_, eigenvectors_->backend());
      } else {
        if (options_.template get<bool>("compute_eigenvalues"))
          eigenvalues_ = std::make_unique<std::vector<XT::Common::complex_t<RealType>>>(
              internal::compute_all_eigenvalues_using_eigen(matrix_.backend()));
        if (options_.template get<bool>("compute_eigenvectors"))
          eigenvectors_ = std::make_unique<EigenDenseMatrix<XT::Common::complex_t<S>>>(
              internal::compute_all_eigenvectors_using_eigen(matrix_.backend()));
      }
#if HAVE_LAPACKE
    } else if (type == "lapack") {
      if (!options_.template get<bool>("compute_eigenvectors"))
        eigenvalues_ = std::make_unique<std::vector<XT::Common::complex_t<RealType>>>(
            internal::compute_eigenvalues_using_lapack(matrix_));
      else {
        eigenvalues_ = std::make_unique<std::vector<XT::Common::complex_t<RealType>>>(N);
        eigenvectors_ = std::make_unique<EigenDenseMatrix<XT::Common::complex_t<S>>>(N, N);
        internal::compute_eigenvalues_and_right_eigenvectors_using_lapack(matrix_, *eigenvalues_, *eigenvectors_);
      }
#endif // HAVE_LAPACKE
    } else
      DUNE_THROW(Common::Exceptions::internal_error,
                 "Given type '" << type << "' is none of EigenSolverOptions<EigenDenseMatrix<S>>::types(), and "
                                           "internal::EigenSolverBase promised to check this!"
                                << "\n\nThese are the available types:\n\n"
                                << EigenSolverOptions<EigenDenseMatrix<S>>::types());
  } // ... compute(...)

  using BaseType::matrix_;
  using BaseType::options_;
  using BaseType::eigenvalues_;
  using BaseType::eigenvectors_;
}; // class EigenSolver<EigenDenseMatrix<...>>


#else // HAVE_EIGEN


template <class S>
class EigenSolverOptions<EigenDenseMatrix<S>>
{
  static_assert(AlwaysFalse<S>::value, "You are missing eigen!");
};


template <class S>
class EigenSolver<EigenDenseMatrix<S>>
{
  static_assert(AlwaysFalse<S>::value, "You are missing eigen!");
};


#endif // HAVE_EIGEN


} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_EIGEN_SOLVER_EIGEN_HH
