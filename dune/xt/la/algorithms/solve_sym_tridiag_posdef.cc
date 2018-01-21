// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)
//   Tobias Leibner  (2017)

#include "config.h"

#include <iostream>
#include <cmath>

#include <dune/common/exceptions.hh>

#include "solve_sym_tridiag_posdef.hh"

#include <dune/xt/common/lapacke.hh>

namespace Dune {
namespace XT {
namespace LA {


void prepare_sym_tridiag_posdef(
    const double* A, size_t dimRange, double* diagonal_elements, double* sub_diagonal_elements, double* anorm)
{
  // get norm of A as ||A||_1 = max_j \sum_i |a_ij|
  *anorm = std::abs(A[0 * dimRange + 0]) + std::abs(A[1 * dimRange + 0]);
  for (size_t jj = 1; jj < dimRange - 1; ++jj) {
    double col_sum = 0.;
    for (size_t ii = jj - 1; ii <= jj + 1; ++ii)
      col_sum += std::abs(A[ii * dimRange + jj]);
    *anorm = std::max(*anorm, col_sum);
  }
  *anorm = std::max(*anorm,
                    std::abs(A[(dimRange - 2) * dimRange + (dimRange - 1)])
                        + std::abs(A[(dimRange - 1) * dimRange + (dimRange - 1)]));

  // get factorization LDL^T
  for (size_t ii = 0; ii < dimRange; ++ii)
    diagonal_elements[ii] = A[ii * dimRange + ii];
  for (size_t ii = 1; ii < dimRange; ++ii)
    sub_diagonal_elements[ii - 1] = A[ii * dimRange + (ii - 1)];
} // void prepare_sym_triadiag_posdef(...)


double solve_sym_tridiag_posdef(
    size_t dimRange, double* diagonal_elements, double* sub_diagonal_elements, double anorm, double* b)
{
  auto info = Common::Lapacke::dpttrf(dimRange, diagonal_elements, sub_diagonal_elements);
  if (info) {
    std::cout << "factorization failed" << std::endl;
    DUNE_THROW(Dune::MathError, "factorization failed");
  }

  // estimate condition of H_k
  double inverse_condition;
  info = Common::Lapacke::dptcon(dimRange, diagonal_elements, sub_diagonal_elements, anorm, &inverse_condition);
  if (info)
    DUNE_THROW(Dune::MathError, "condition estimation failed");
  if (inverse_condition < 1e-20) {
    std::cout << "inverse_condition_failed" << std::endl;
    DUNE_THROW(Dune::MathError, "matrix singular");
  }

  // solve system
  info = Common::Lapacke::dpttrs(
      Common::Lapacke::row_major(), dimRange, 1, diagonal_elements, sub_diagonal_elements, b, 1);
  if (info) {
    std::cout << "solving system failed" << std::endl;
    DUNE_THROW(Dune::MathError, "solving system failed");
  }
  return inverse_condition;
}


} // namespace LA
} // namespace GDT
} // namespace Dune
