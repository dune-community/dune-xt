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

namespace Dune {
namespace XT {
namespace LA {


/**
 * \brief Extracts input needed for solve_sym_tridiag_posdef from matrix A
 */
void prepare_sym_tridiag_posdef(
    const double* A, size_t dimRange, double* diagonal_elements, double* sub_diagonal_elements, double* anorm);

/**  \brief Solves linear equation for tridiagonal symmetric positive definite matrix.
  * \param[in] dimRange number of rows and columns of matrix
  * \param[in] diagonal_elements array containing the diagonal elements of the matrix (length dimRange)
  * \param[in] sub_diagonal_elements array containing the sub-diagonal elements of the matrix (length dimRange-1)
  * \param[in] anorm operator norm of the matrix (1-norm)
  * \param[in/out] b array containing the rhs (length dimRange). Is overwritten by the solution of the equation.
  * \returns estimate of inverse of the condition of the matrix
  * \attention This function depends on LAPACKE. If LAPACKE is not found an error is thrown.
  */
double solve_sym_tridiag_posdef(
    size_t dimRange, double* diagonal_elements, double* sub_diagonal_elements, double anorm, double* b);


} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_ALGORITHMS_SOLVE_SYM_TRIDIAG_POSDEF_HH
