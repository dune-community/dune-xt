// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017 - 2018)

#ifndef DUNE_XT_LA_ALGORITHMS_SOLVE_SYM_TRIDIAG_POSDEF_HH
#define DUNE_XT_LA_ALGORITHMS_SOLVE_SYM_TRIDIAG_POSDEF_HH

#include <cstddef>

#include <dune/common/exceptions.hh>

#include <dune/xt/common/lapacke.hh>
#include <dune/xt/common/matrix.hh>
#include <dune/xt/common/vector.hh>

#include <dune/xt/la/algorithms/cholesky.hh>

namespace Dune {
namespace XT {
namespace LA {


/**  \brief Solves linear equation for tridiagonal symmetric positive definite matrix.
  * \param[in] dimRange number of rows and columns of matrix
  * \param[in] diag array containing the diagonal elements of the matrix (length dimRange)
  * \param[in] sub_diag array containing the sub-diagonal elements of the matrix (length dimRange-1)
  * \param[in] anorm operator norm of the matrix (1-norm)
  * \param[in/out] b array containing the rhs (length dimRange). Is overwritten by the solution of the equation.
  * \returns estimate of inverse of the condition of the matrix
  * \attention This function depends on LAPACKE. If LAPACKE is not found an error is thrown.
  */
template <class VectorType, class SecondVectorType, class RhsVectorType>
std::enable_if_t<Common::is_vector<VectorType>::value && Common::is_vector<SecondVectorType>::value
                     && Common::is_vector<RhsVectorType>::value,
                 void>
solve_sym_tridiag_posdef(VectorType& diag, SecondVectorType& sub_diag, RhsVectorType& b)
{
  tridiagonal_cholesky(diag, subdiag);
  solve_tridiagonal_cholesky_factorized(diag, subdiag, b);
}

template <class MatrixType, class VectorType, class RhsVectorType>
std::enable_if_t<Common::is_matrix<MatrixType>::value && Common::is_vector<VectorType>::value
                     && Common::is_vector<RhsVectorType>::value,
                 void>
solve_sym_tridiag_posdef(const MatrixType& A, VectorType& x, const RhsVectorType& y)
{
  typedef Common::MatrixAbstraction<MatrixType> Mat;
  typedef Common::VectorAbstraction<VectorType> Vec;
  typedef Common::VectorAbstraction<RhsVectorType> RhsVec;
  typedef typename Mat::ScalarType ScalarType;
  const size_t num_rows = Mat::rows(A);
  std::vector<ScalarType> diag(num_rows, 0.);
  std::vector<ScalarType> sub_diag(num_rows - 1, 0.);
  for (size_t rr = 0; rr < num_rows; ++rr)
    diag[rr] = Mat::get_entry(A, rr, rr);
  for (size_t rr = 0; rr < num_rows - 1; ++rr)
    sub_diag[rr] = Mat::get_entry(A, rr + 1, rr);
  // copy y to x
  x = Common::convert_to<VectorType>(y);
  solve_sym_tridiag_posdef(diag, sub_diag, x);
}


} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_ALGORITHMS_SOLVE_SYM_TRIDIAG_POSDEF_HH
