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

#ifndef DUNE_XT_LA_EIGEN_SOLVER_FMATRIX_HH
#define DUNE_XT_LA_EIGEN_SOLVER_FMATRIX_HH

#include <algorithm>
#include <functional>

#include <dune/common/typetraits.hh>

#include <dune/xt/common/fmatrix.hh>
#include <dune/xt/common/matrix.hh>
#include <dune/xt/common/numeric_cast.hh>
#include <dune/xt/common/vector.hh>
#include <dune/xt/la/container/conversion.hh>
#include <dune/xt/la/container/eigen/dense.hh>
#include <dune/xt/la/solver.hh>
#include <dune/xt/la/eigen-solver.hh>

#include "internal/base.hh"
#include "internal/eigen.hh"
#include "internal/lapacke.hh"
#include "internal/numpy.hh"
#include "internal/shifted-qr.hh"

namespace Dune {
namespace XT {
namespace LA {


template <class K, int SIZE>
class EigenSolverOptions<Dune::FieldMatrix<K, SIZE, SIZE>>
{
public:
  static std::vector<std::string> types()
  {
    std::vector<std::string> tps;
    if (Common::Lapacke::available())
      tps.push_back("lapack");
#if HAVE_EIGEN
    tps.push_back("eigen");
#endif
    if (internal::numpy_eigensolver_available())
      tps.push_back("numpy");
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
}; // class EigenSolverOptions<Dune::FieldMatrix<K, SIZE, SIZE>>


template <class K, int SIZE>
class EigenSolverOptions<Dune::XT::Common::FieldMatrix<K, SIZE, SIZE>>
    : public EigenSolverOptions<Dune::FieldMatrix<K, SIZE, SIZE>>
{
};


template <class K, int SIZE>
class EigenSolver<Dune::FieldMatrix<K, SIZE, SIZE>>
    : public internal::EigenSolverBase<Dune::FieldMatrix<K, SIZE, SIZE>,
                                       K,
                                       Dune::FieldMatrix<XT::Common::real_t<K>, SIZE, SIZE>,
                                       Dune::FieldMatrix<XT::Common::complex_t<K>, SIZE, SIZE>>
{
  using BaseType = internal::EigenSolverBase<Dune::FieldMatrix<K, SIZE, SIZE>,
                                             K,
                                             Dune::FieldMatrix<XT::Common::real_t<K>, SIZE, SIZE>,
                                             Dune::FieldMatrix<XT::Common::complex_t<K>, SIZE, SIZE>>;

public:
  using typename BaseType::MatrixType;

  template <class... Args>
  explicit EigenSolver(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }

protected:
  void compute() const override final
  {
    const auto type = options_->template get<std::string>("type");
#if HAVE_LAPACKE
    if (type == "lapack") {
      if (!options_->template get<bool>("compute_eigenvectors"))
        eigenvalues_ = std::make_unique<std::vector<XT::Common::complex_t<K>>>(
            internal::compute_eigenvalues_using_lapack(matrix_));
      else {
        eigenvalues_ = std::make_unique<std::vector<XT::Common::complex_t<K>>>(SIZE);
        eigenvectors_ = std::make_unique<Dune::FieldMatrix<XT::Common::complex_t<K>, SIZE, SIZE>>();
        internal::compute_eigenvalues_and_right_eigenvectors_using_lapack(matrix_, *eigenvalues_, *eigenvectors_);
      }
    } else
#endif // HAVE_LAPACKE
#if HAVE_EIGEN
        if (type == "eigen") {
      if (options_->template get<bool>("compute_eigenvalues") && options_->template get<bool>("compute_eigenvectors")) {
        eigenvalues_ = std::make_unique<std::vector<XT::Common::complex_t<K>>>(SIZE);
        EigenDenseMatrix<K> tmp_matrix(matrix_);
        EigenDenseMatrix<Common::complex_t<K>> tmp_eigenvectors(matrix_);
        internal::compute_eigenvalues_and_right_eigenvectors_using_eigen(
            tmp_matrix.backend(), *eigenvalues_, tmp_eigenvectors.backend());
        eigenvectors_ = std::make_unique<Dune::FieldMatrix<XT::Common::complex_t<K>, SIZE, SIZE>>(
            convert_to<Dune::FieldMatrix<XT::Common::complex_t<K>, SIZE, SIZE>>(tmp_eigenvectors));
      } else {
        if (options_->template get<bool>("compute_eigenvalues"))
          eigenvalues_ = std::make_unique<std::vector<XT::Common::complex_t<K>>>(
              internal::compute_eigenvalues_using_eigen(EigenDenseMatrix<K>(matrix_).backend()));
        if (options_->template get<bool>("compute_eigenvectors")) {
          eigenvectors_ = std::make_unique<Dune::FieldMatrix<XT::Common::complex_t<K>, SIZE, SIZE>>(
              convert_to<Dune::FieldMatrix<XT::Common::complex_t<K>, SIZE, SIZE>>(
                  EigenDenseMatrix<XT::Common::complex_t<K>>(
                      internal::compute_right_eigenvectors_using_eigen(EigenDenseMatrix<K>(matrix_).backend()))));
        }
      }
    } else
#endif // HAVE_EIGEN
        if (type == "numpy") {
      if (options_->template get<bool>("compute_eigenvalues") || options_->template get<bool>("compute_eigenvectors")) {
        eigenvalues_ = std::make_unique<std::vector<XT::Common::complex_t<K>>>(SIZE);
        eigenvectors_ = std::make_unique<Dune::FieldMatrix<XT::Common::complex_t<K>, SIZE, SIZE>>();
        internal::compute_eigenvalues_and_right_eigenvectors_of_a_fieldmatrix_using_numpy(
            matrix_, *eigenvalues_, *eigenvectors_);
      }
    } else if (type == "shifted_qr") {
      if (options_->template get<bool>("compute_eigenvalues") || options_->template get<bool>("compute_eigenvectors")) {
        eigenvalues_ = std::make_unique<std::vector<XT::Common::complex_t<K>>>(SIZE);
        eigenvectors_ = std::make_unique<Dune::FieldMatrix<XT::Common::complex_t<K>, SIZE, SIZE>>();
        std::vector<XT::Common::real_t<K>> real_eigenvalues(SIZE);
        auto real_eigenvectors = std::make_unique<Dune::FieldMatrix<XT::Common::real_t<K>, SIZE, SIZE>>();
        internal::compute_real_eigenvalues_and_real_right_eigenvectors_using_qr(
            matrix_, real_eigenvalues, *real_eigenvectors);
        for (size_t ii = 0; ii < SIZE; ++ii) {
          (*eigenvalues_)[ii] = real_eigenvalues[ii];
          for (size_t jj = 0; jj < SIZE; ++jj)
            (*eigenvectors_)[ii][jj] = (*real_eigenvectors)[ii][jj];
        }
      }
    } else
      DUNE_THROW(Common::Exceptions::internal_error,
                 "Given type '" << type << "' is none of EigenSolverOptions<Dune::FieldMatrix<K, ROWS, "
                                           "COLS>>::types(), and  internal::EigenSolverBase promised to check this!"
                                << "\n\nThese are the available types:\n\n"
                                << EigenSolverOptions<MatrixType>::types());
  } //... compute(...)

  using BaseType::matrix_;
  using BaseType::options_;
  using BaseType::eigenvalues_;
  using BaseType::eigenvectors_;
}; // class EigenSolver<FieldMatrix<...>>


template <class K, int SIZE>
class EigenSolver<Dune::XT::Common::FieldMatrix<K, SIZE, SIZE>> : public EigenSolver<Dune::FieldMatrix<K, SIZE, SIZE>>
{
public:
  template <class... Args>
  EigenSolver(Args&&... args)
    : EigenSolver<Dune::FieldMatrix<K, SIZE, SIZE>>(std::forward<Args>(args)...)
  {
  }
};


} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_EIGEN_SOLVER_FMATRIX_HH
