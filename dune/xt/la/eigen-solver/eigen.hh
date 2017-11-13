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
class EigenDenseMatrixEigenSolverTraits
{
public:
  typedef EigenDenseMatrix<S> MatrixType;
  typedef typename MatrixType::RealType RealType;
  typedef typename std::complex<RealType> ComplexType;
  typedef typename XT::LA::Container<RealType, MatrixType::vector_type>::VectorType RealVectorType;
  typedef typename XT::LA::Container<ComplexType, MatrixType::vector_type>::VectorType ComplexVectorType;
  typedef EigenDenseMatrix<RealType> RealMatrixType;
  typedef EigenDenseMatrix<ComplexType> ComplexMatrixType;
  typedef EigenSolver<MatrixType> derived_type;
};

template <class S>
class EigenSolver<EigenDenseMatrix<S>> : public EigenSolverBase<EigenDenseMatrixEigenSolverTraits<S>>
{
  typedef EigenSolverBase<EigenDenseMatrixEigenSolverTraits<S>> BaseType;

public:
  using typename BaseType::MatrixType;
  using typename BaseType::ComplexType;
  using typename BaseType::ComplexVectorType;

  EigenSolver(const MatrixType& matrix)
    : BaseType(matrix)
  {
  }

  static std::vector<std::string> types()
  {
    return
    {
      "eigen",
#if HAVE_LAPACKE
          "lapacke",
#endif
#if 0
      ,
      "qrhouseholder"
#endif
    };
  }

  static Common::Configuration options(const std::string type = "")
  {
    const std::string tp = !type.empty() ? type : types()[0];
    internal::SolverUtils::check_given(tp, types());
    return {{"type", tp},
            {"check_for_inf_nan", "1"},
            {"check_evs_are_real", "0"},
            {"check_evs_are_positive", "0"},
            {"check_eigenvectors_are_real", "0"}};
  }

  virtual void get_eigenvalues(std::vector<ComplexType>& evs, const std::string& type) const override final
  {
    if (type == "eigen")
      evs = internal::compute_all_eigenvalues_using_eigen(matrix_.backend());
#if HAVE_LAPACKE
    else if (type == "lapacke")
      evs = internal::compute_all_eigenvalues_using_lapacke(matrix_);
#endif
    else
      DUNE_THROW(Common::Exceptions::internal_error,
                 "Given type '" << type << "' is not supported, although it was reported by types()!");
  }

  virtual void get_eigenvectors(std::vector<ComplexVectorType>& evs, const std::string& type) const override final
  {
    if (type == "eigen")
      evs = internal::compute_all_eigenvectors_using_eigen(matrix_.backend());
#if HAVE_LAPACKE
    else if (type == "lapacke")
      internal::compute_all_eigenvectors_using_lapacke(matrix_, evs, evs[0][0]);
#endif
    else
      DUNE_THROW(Common::Exceptions::internal_error,
                 "Given type '" << type << "' is not supported, although it was reported by types()!");
  }

protected:
  using BaseType::matrix_;
}; // class EigenSolver<EigenDenseMatrix<S>>


#else // HAVE_EIGEN


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
