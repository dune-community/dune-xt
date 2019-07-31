// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#ifndef DUNE_XT_LA_GENERALIZED_EIGEN_SOLVER_DEFAULT_HH
#define DUNE_XT_LA_GENERALIZED_EIGEN_SOLVER_DEFAULT_HH

#include <vector>

#include <dune/xt/common/matrix.hh>
#include <dune/xt/la/container/conversion.hh>
#include <dune/xt/la/eigen-solver.hh>
#include <dune/xt/la/exceptions.hh>

#include "internal/base.hh"
#include "internal/lapacke.hh"

namespace Dune {
namespace XT {
namespace LA {


template <class MatrixType>
class GeneralizedEigenSolverOptions<MatrixType, true>
{
public:
  static std::vector<std::string> types()
  {
    std::vector<std::string> tps;
    if (Common::Lapacke::available())
      tps.push_back("lapack");
    DUNE_THROW_IF(tps.empty(),
                  Exceptions::generalized_eigen_solver_failed,
                  "No backend available for generalized eigenvalue problems!");
    return tps;
  }

  static Common::Configuration options(const std::string type = "")
  {
    const std::string actual_type = type.empty() ? types()[0] : type;
    internal::ensure_generalized_eigen_solver_type(actual_type, types());
    Common::Configuration opts = internal::default_generalized_eigen_solver_options();
    opts["type"] = actual_type;
    return opts;
  }
}; // class GeneralizedEigenSolverOptions<>


template <class MatrixImp>
class GeneralizedEigenSolver<MatrixImp, true>
  : public internal::GeneralizedEigenSolverBase<
        MatrixImp,
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
  using BaseType = internal::GeneralizedEigenSolverBase<
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
  explicit GeneralizedEigenSolver(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {}

protected:
  void compute() const override final
  {
    const auto type = options_->template get<std::string>("type");
    if (type == "lapack") {
      DUNE_THROW_IF(!Common::Lapacke::available(),
                    Exceptions::generalized_eigen_solver_failed_bc_it_was_not_set_up_correctly,
                    "Lapacke backend not available!");
      if (!options_->template get<bool>("compute_eigenvectors"))
        eigenvalues_ = std::make_unique<std::vector<ComplexType>>(
            internal::compute_generalized_eigenvalues_using_lapack(lhs_matrix_, rhs_matrix_));
      else {
        DUNE_THROW(Exceptions::generalized_eigen_solver_failed,
                   "Eigenvectors of generalized eigenvalue problems are not yet implemented using lapacke!");
      }
    } else
      DUNE_THROW(Common::Exceptions::internal_error,
                 "Given type '" << type << "' is none of GeneralizedEigenSolverOptions<...>::types(),"
                                << " and internal::GeneralizedEigenSolverBase promised to check this!"
                                << "\n\nThese are the available types:\n\n"
                                << GeneralizedEigenSolverOptions<MatrixType>::types());
  } //... compute(...)

  using BaseType::eigenvalues_;
  using BaseType::eigenvectors_;
  using BaseType::lhs_matrix_;
  using BaseType::options_;
  using BaseType::rhs_matrix_;
}; // class GeneralizedEigenSolver<MatrixType, true>


} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_GENERALIZED_EIGEN_SOLVER_DEFAULT_HH
