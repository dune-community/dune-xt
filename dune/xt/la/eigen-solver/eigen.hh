// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2017 - 2018)
//   Tobias Leibner  (2017 - 2018)

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
#include "internal/shifted-qr.hh"

namespace Dune {
namespace XT {
namespace LA {

#if HAVE_EIGEN


template <class S>
class EigenSolverOptions<EigenDenseMatrix<S>, true>
{
public:
  static std::vector<std::string> types()
  {
    std::vector<std::string> tps = {"eigen"};
    if (Common::Lapacke::available())
      tps.push_back("lapack");
    tps.push_back("shifted_qr");
    return tps;
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
class EigenSolver<EigenDenseMatrix<S>, true>
  : public internal::EigenSolverBase<EigenDenseMatrix<S>,
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
  {}

protected:
  void compute() const override final
  {
    const auto type = options_->template get<std::string>("type", EigenSolverOptions<EigenDenseMatrix<S>>::types()[0]);
    const size_t N = matrix_.rows();
    const bool compute_eigenvalues = options_->template get<bool>("compute_eigenvalues", true);
    const bool compute_eigenvectors = options_->template get<bool>("compute_eigenvectors", false);
    if (type == "eigen") {
      if (compute_eigenvalues && compute_eigenvectors) {
        eigenvalues_ = std::make_unique<std::vector<XT::Common::complex_t<RealType>>>(N);
        eigenvectors_ = std::make_unique<EigenDenseMatrix<XT::Common::complex_t<S>>>(N, N);
        internal::compute_eigenvalues_and_right_eigenvectors_using_eigen(
            matrix_.backend(), *eigenvalues_, eigenvectors_->backend());
      } else if (compute_eigenvalues)
        eigenvalues_ = std::make_unique<std::vector<XT::Common::complex_t<RealType>>>(
            internal::compute_eigenvalues_using_eigen(matrix_.backend()));
      else if (compute_eigenvectors)
        eigenvectors_ = std::make_unique<EigenDenseMatrix<XT::Common::complex_t<S>>>(
            internal::compute_right_eigenvectors_using_eigen(matrix_.backend()));
#if HAVE_LAPACKE || HAVE_MKL
    } else if (type == "lapack") {
      if (!compute_eigenvectors)
        eigenvalues_ = std::make_unique<std::vector<XT::Common::complex_t<RealType>>>(
            internal::compute_eigenvalues_using_lapack(matrix_));
      else {
        eigenvalues_ = std::make_unique<std::vector<XT::Common::complex_t<RealType>>>(N);
        eigenvectors_ = std::make_unique<EigenDenseMatrix<XT::Common::complex_t<S>>>(N, N);
        internal::compute_eigenvalues_and_right_eigenvectors_using_lapack(matrix_, *eigenvalues_, *eigenvectors_);
      }
#endif // HAVE_LAPACKE || HAVE_MKL
    } else if (type == "shifted_qr") {
      if (!compute_eigenvectors) {
        eigenvalues_ = std::make_unique<std::vector<XT::Common::complex_t<RealType>>>(N);
        std::vector<XT::Common::real_t<RealType>> real_eigenvalues(N);
        real_eigenvalues = internal::compute_eigenvalues_using_qr(matrix_);
        for (size_t ii = 0; ii < N; ++ii)
          (*eigenvalues_)[ii] = real_eigenvalues[ii];
      } else {
        eigenvalues_ = std::make_unique<std::vector<XT::Common::complex_t<RealType>>>(N);
        eigenvectors_ = std::make_unique<EigenDenseMatrix<XT::Common::complex_t<S>>>(N, N);
        std::vector<XT::Common::real_t<RealType>> real_eigenvalues(N);
        auto real_eigenvectors = std::make_unique<Dune::DynamicMatrix<XT::Common::real_t<RealType>>>(N, N);
        internal::compute_real_eigenvalues_and_real_right_eigenvectors_using_qr(
            matrix_, real_eigenvalues, *real_eigenvectors);
        for (size_t ii = 0; ii < N; ++ii) {
          (*eigenvalues_)[ii] = real_eigenvalues[ii];
          for (size_t jj = 0; jj < N; ++jj)
            (*eigenvectors_).set_entry(ii, jj, (*real_eigenvectors)[ii][jj]);
        }
      }
    } else
      DUNE_THROW(Common::Exceptions::internal_error,
                 "Given type '" << type
                                << "' is none of EigenSolverOptions<EigenDenseMatrix<S>>::types(), and "
                                   "internal::EigenSolverBase promised to check this!"
                                << "\n\nThese are the available types:\n\n"
                                << EigenSolverOptions<EigenDenseMatrix<S>>::types());
  } // ... compute(...)

  using BaseType::eigenvalues_;
  using BaseType::eigenvectors_;
  using BaseType::matrix_;
  using BaseType::options_;
}; // class EigenSolver<EigenDenseMatrix<...>>


#else // HAVE_EIGEN


template <class S>
class EigenSolverOptions<EigenDenseMatrix<S>, true>
{
  static_assert(AlwaysFalse<S>::value, "You are missing eigen!");
};


template <class S>
class EigenSolver<EigenDenseMatrix<S>, true>
{
  static_assert(AlwaysFalse<S>::value, "You are missing eigen!");
};


#endif // HAVE_EIGEN


} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_EIGEN_SOLVER_EIGEN_HH
