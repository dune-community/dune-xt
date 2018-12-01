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

#ifndef DUNE_XT_LA_EIGEN_SOLVER_DEFAULT_HH
#define DUNE_XT_LA_EIGEN_SOLVER_DEFAULT_HH

#include "config.h"

#include <vector>

#include <dune/xt/common/matrix.hh>
#include <dune/xt/la/container/conversion.hh>
#include <dune/xt/la/eigen-solver.hh>

#include "internal/base.hh"
#include "internal/shifted-qr.hh"
#include "internal/lapacke.hh"

namespace Dune {
namespace XT {
namespace LA {


template <class MatrixType>
class EigenSolverOptions<MatrixType, true>
{
public:
  static std::vector<std::string> types()
  {
    std::vector<std::string> tps;
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
}; // class EigenSolverOptions<>


template <class MatrixImp>
class EigenSolver<MatrixImp, true>
  : public internal::EigenSolverBase<MatrixImp,
                                     typename Common::MatrixAbstraction<MatrixImp>::ScalarType,
                                     typename Common::MatrixAbstraction<MatrixImp>::template MatrixTypeTemplate<
                                         Common::MatrixAbstraction<MatrixImp>::static_rows,
                                         Common::MatrixAbstraction<MatrixImp>::static_cols,
                                         typename Common::MatrixAbstraction<MatrixImp>::RealType>,
                                     typename Common::MatrixAbstraction<MatrixImp>::template MatrixTypeTemplate<
                                         Common::MatrixAbstraction<MatrixImp>::static_rows,
                                         Common::MatrixAbstraction<MatrixImp>::static_cols,
                                         std::complex<typename Common::MatrixAbstraction<MatrixImp>::RealType>>>
{
  using M = Common::MatrixAbstraction<MatrixImp>;
  using BaseType = internal::EigenSolverBase<
      MatrixImp,
      typename M::ScalarType,
      typename M::template MatrixTypeTemplate<M::static_rows, M::static_cols, typename M::RealType>,
      typename M::template MatrixTypeTemplate<M::static_rows, M::static_cols, std::complex<typename M::RealType>>>;

public:
  using typename BaseType::ComplexMatrixType;
  using typename BaseType::ComplexType;
  using typename BaseType::MatrixType;
  using typename BaseType::RealMatrixType;
  using typename BaseType::RealType;
  using RealM = Common::MatrixAbstraction<RealMatrixType>;
  using ComplexM = Common::MatrixAbstraction<ComplexMatrixType>;

  template <class... Args>
  explicit EigenSolver(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {}

protected:
  void compute() const override final
  {
    const auto type = options_->template get<std::string>("type");
    const auto rows = M::rows(matrix_);
    const auto cols = M::rows(matrix_);
#if HAVE_LAPACKE || HAVE_MKL
    if (type == "lapack") {
      if (!options_->template get<bool>("compute_eigenvectors"))
        eigenvalues_ = std::make_unique<std::vector<ComplexType>>(internal::compute_eigenvalues_using_lapack(matrix_));
      else {
        eigenvalues_ = std::make_unique<std::vector<ComplexType>>(rows);
        eigenvectors_ = ComplexM::make_unique(rows, cols);
        internal::compute_eigenvalues_and_right_eigenvectors_using_lapack(matrix_, *eigenvalues_, *eigenvectors_);
      }
    } else
#endif // HAVE_LAPACKE || HAVE_MKL
        if (type == "shifted_qr") {
      if (options_->template get<bool>("compute_eigenvalues") || options_->template get<bool>("compute_eigenvectors")) {
        eigenvalues_ = std::make_unique<std::vector<ComplexType>>(rows);
        eigenvectors_ = ComplexM::make_unique(rows, cols);
        std::vector<RealType> real_eigenvalues(rows);
        auto real_eigenvectors = RealM::make_unique(rows, cols);
        internal::compute_real_eigenvalues_and_real_right_eigenvectors_using_qr(
            matrix_, real_eigenvalues, *real_eigenvectors);
        for (size_t ii = 0; ii < rows; ++ii) {
          (*eigenvalues_)[ii] = real_eigenvalues[ii];
          for (size_t jj = 0; jj < cols; ++jj)
            ComplexM::set_entry(*eigenvectors_, ii, jj, RealM::get_entry(*real_eigenvectors, ii, jj));
        }
      }
    } else
      DUNE_THROW(Common::Exceptions::internal_error,
                 "Given type '" << type << "' is none of EigenSolverOptions<...>::types(),"
                                << " and internal::EigenSolverBase promised to check this!"
                                << "\n\nThese are the available types:\n\n"
                                << EigenSolverOptions<MatrixType>::types());
  } //... compute(...)

  using BaseType::eigenvalues_;
  using BaseType::eigenvectors_;
  using BaseType::matrix_;
  using BaseType::options_;
}; // class EigenSolver<MatrixType, true>


} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_EIGEN_SOLVER_DEFAULT_HH
