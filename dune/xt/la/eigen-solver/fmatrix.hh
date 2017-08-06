// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2017 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_XT_LA_EIGEN_SOLVER_FMATRIX_HH
#define DUNE_XT_LA_EIGEN_SOLVER_FMATRIX_HH

#include <algorithm>
#include <functional>

#include <dune/xt/la/container/eigen/dense.hh>
#include <dune/xt/la/solver.hh>

#include "../eigen-solver.hh"
#include "internal/lapacke.hh"

namespace Dune {
namespace XT {
namespace LA {


template <class S, int dimRange>
class FieldMatrixEigenSolverTraits
{
public:
  typedef Dune::FieldMatrix<S, dimRange, dimRange> MatrixType;
  typedef typename Dune::FieldTraits<S>::real_type RealType;
  typedef typename std::complex<RealType> ComplexType;
  typedef typename Dune::FieldVector<RealType, dimRange> RealVectorType;
  typedef typename Dune::FieldVector<ComplexType, dimRange> ComplexVectorType;
  typedef typename Dune::FieldMatrix<RealType, dimRange, dimRange> RealMatrixType;
  typedef typename Dune::FieldMatrix<ComplexType, dimRange, dimRange> ComplexMatrixType;
  typedef EigenSolver<MatrixType> derived_type;
};

template <class S, int dimRange>
class EigenSolver<Dune::FieldMatrix<S, dimRange, dimRange>>
    : public EigenSolverBase<FieldMatrixEigenSolverTraits<S, dimRange>>
{
  typedef EigenSolverBase<FieldMatrixEigenSolverTraits<S, dimRange>> BaseType;

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
#if HAVE_LAPACKE
      "lapacke",
#endif
#if HAVE_EIGEN
          "eigen",
#endif
          "qrhouseholder"
    };
  }

  static Common::Configuration options(const std::string type = "")
  {
    const std::string tp = !type.empty() ? type : types()[0];
    internal::SolverUtils::check_given(tp, types());
    Common::Configuration default_options(
        {"type", "check_for_inf_nan", "check_evs_are_real", "check_evs_are_positive", "check_eigenvectors_are_real"},
        {tp, "1", "0", "0", "0"});
    return default_options;
  }

  virtual void get_eigenvalues(std::vector<ComplexType>& evs, const std::string& type) const override final
  {
#if HAVE_LAPACKE
    if (type == "lapacke")
      evs = internal::compute_all_eigenvalues_using_lapacke(matrix_);
#endif
#if HAVE_EIGEN
    else if (type == "eigen") {
      XT::LA::EigenDenseMatrix<S> eigenmatrix(dimRange, dimRange);
      for (size_t rr = 0; rr < dimRange; ++rr)
        for (size_t cc = 0; cc < dimRange; ++cc)
          eigenmatrix.set_entry(rr, cc, matrix_[rr][cc]);
      evs = internal::compute_all_eigenvalues_using_eigen(eigenmatrix.backend());
    }
#endif
    else
      DUNE_THROW(Common::Exceptions::internal_error,
                 "Given type '" << type << "' is not supported, although it was reported by types()!");
  } // ... get_eigenvalues(...)

  virtual void get_eigenvectors(std::vector<ComplexVectorType>& evs, const std::string& type) const override final
  {
    evs.resize(dimRange);
#if HAVE_LAPACKE
    if (type == "lapacke")
      internal::compute_all_eigenvectors_using_lapacke(matrix_, evs, evs[0][0]);
#endif
#if HAVE_EIGEN
    else if (type == "eigen") {
      XT::LA::EigenDenseMatrix<S> eigenmatrix(dimRange, dimRange);
      for (size_t rr = 0; rr < dimRange; ++rr)
        for (size_t cc = 0; cc < dimRange; ++cc)
          eigenmatrix.set_entry(rr, cc, matrix_[rr][cc]);
      auto eigen_eigvecs = internal::compute_all_eigenvectors_using_eigen(eigenmatrix.backend());
      for (size_t rr = 0; rr < dimRange; ++rr)
        for (size_t cc = 0; cc < dimRange; ++cc)
          evs[rr][cc] = eigen_eigvecs[rr].get_entry(cc);
    }
#endif
    else
      DUNE_THROW(Common::Exceptions::internal_error,
                 "Given type '" << type << "' is not supported, although it was reported by types()!");
  } // ... get_eigenvectors(...)

private:
  using BaseType::matrix_;
}; // class EigenSolver<FieldMatrix<S>>


} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_EIGEN_SOLVER_FMATRIX_HH
