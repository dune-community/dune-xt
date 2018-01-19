// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_XT_LA_ALGORITHMS_QR_HH
#define DUNE_XT_LA_ALGORITHMS_QR_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/math.hh>

#include <dune/xt/common/float_cmp.hh>

#include <dune/xt/la/container.hh>

namespace Dune {
namespace XT {
namespace LA {


// Calculates A = H * A(row_begin:row_end,col_begin:col_end) where H = I-tau*w*w^T and w = v[row_begin:row_end]
template <class FieldType, int rows, int cols>
void multiply_householder_from_left(FieldMatrix<FieldType, rows, cols>& A,
                                    const FieldType& tau,
                                    const FieldVector<FieldType, rows>& v,
                                    const size_t row_begin,
                                    const size_t row_end,
                                    const size_t col_begin,
                                    const size_t col_end)
{
  // calculate w^T A first
  FieldVector<FieldType, cols> wT_A(0.);
  for (size_t rr = row_begin; rr < row_end; ++rr)
    for (size_t cc = col_begin; cc < col_end; ++cc)
      wT_A[cc] += v[rr] * A[rr][cc];
  for (size_t rr = row_begin; rr < row_end; ++rr)
    for (size_t cc = col_begin; cc < col_end; ++cc)
      A[rr][cc] -= tau * v[rr] * wT_A[cc];
}

// Calculates A = A * H.
// \see multiply_householder_from_left
template <class FieldType, int rows, int cols>
void multiply_householder_from_right(FieldMatrix<FieldType, rows, cols>& A,
                                     const FieldType& tau,
                                     const FieldVector<FieldType, cols>& v,
                                     const size_t row_begin,
                                     const size_t row_end,
                                     const size_t col_begin,
                                     const size_t col_end)
{
  // calculate A w first
  FieldVector<FieldType, rows> Aw(0.);
  for (size_t rr = row_begin; rr < row_end; ++rr)
    for (size_t cc = col_begin; cc < col_end; ++cc)
      Aw[rr] += A[rr][cc] * v[cc];
  for (size_t rr = row_begin; rr < row_end; ++rr)
    for (size_t cc = col_begin; cc < col_end; ++cc)
      A[rr][cc] -= tau * Aw[rr] * v[cc];
}

/** \brief This is a simple QR scheme using Householder reflections and column pivoting.
  * The householder matrix applied in each step is H = I - 2 v v^T, where v = u/||u|| and u = x - s ||x|| e_1,
  * s = +-1 has the opposite sign of u_1 and x is the current column of A. The matrix H is rewritten as
  * H = I - tau w w^T, where w=u/u_1 and tau = -s u_1/||x||.
  * Calculates AP = QR, where P is a permutation matrix reordering the columns of A.
  * The matrix A is overwritten with the QR decomposition in compressed form, i.e. after completion of this function
  * the upper triangular part of A contains R and each column jj of the strictly lower triangular part of A contains
  * w_jj (except for the first element of w_jj, which always is 1). The vector tau contains the tau_jj used for the
  * householder matrices. From this information, Q can be reconstructed if necessary.
  * \see https://en.wikipedia.org/wiki/QR_decomposition#Using_Householder_reflections.
  * \see http://www.cs.cornell.edu/~bindel/class/cs6210-f09/lec18.pdf
  */
template <class FieldType, int rows, int cols>
void qr_decomposition(FieldMatrix<FieldType, rows, cols>& A,
                      FieldVector<FieldType, cols>& tau,
                      FieldVector<size_t, cols>& permutations)
{
  std::fill(tau.begin(), tau.end(), 0.);
  for (size_t ii = 0; ii < cols; ++ii)
    permutations[ii] = ii;

  // compute (squared) column norms
  FieldVector<FieldType, cols> col_norms(0.);
  for (size_t rr = 0; rr < rows; ++rr)
    for (size_t cc = 0; cc < cols; ++cc)
      col_norms[cc] += std::pow(A[rr][cc], 2);

  FieldVector<FieldType, rows> w(0.);

  for (size_t jj = 0; jj < cols; ++jj) {

    // Pivoting
    // swap column jj and column with greatest norm
    auto max_it = std::max_element(col_norms.begin() + jj, col_norms.end());
    size_t max_index = std::distance(col_norms.begin(), max_it);
    if (max_index != jj) {
      for (size_t rr = 0; rr < rows; ++rr)
        std::swap(A[rr][jj], A[rr][max_index]);
      std::swap(col_norms[jj], col_norms[max_index]);
      std::swap(permutations[jj], permutations[max_index]);
    }

    // Matrix update
    // Reduction by householder matrix
    FieldType normx(0);
    for (size_t rr = jj; rr < rows; ++rr)
      normx += std::pow(A[rr][jj], 2);
    normx = std::sqrt(normx);

    if (XT::Common::FloatCmp::ne(normx, 0.)) {
      const auto s = -sign(A[jj][jj]);
      const FieldType u1 = A[jj][jj] - s * normx;
      w[jj] = 1.;
      for (size_t rr = jj + 1; rr < rows; ++rr) {
        w[rr] = A[rr][jj] / u1;
        A[rr][jj] = w[rr];
      }
      A[jj][jj] = s * normx;
      tau[jj] = -s * u1 / normx;
      // calculate A = H A
      multiply_householder_from_left(A, tau[jj], w, jj, rows, jj + 1, cols);
    } // if (normx != 0)

    // Norm downdate
    for (size_t rr = jj + 1; rr < rows; ++rr)
      col_norms[rr] -= std::pow(A[jj][rr], 2);

  } // jj
} // void QR(...)

/**
   * \brief Same QR scheme as above, but Q is explicitly calculated
   * \see QR
   */
template <class FieldType, int rows, int cols>
void qr_decomposition(FieldMatrix<FieldType, rows, cols>& A,
                      FieldVector<size_t, cols>& permutations,
                      FieldMatrix<FieldType, rows, rows>& Q)
{
  std::fill(Q.begin(), Q.end(), 0.);
  for (size_t ii = 0; ii < rows; ++ii)
    Q[ii][ii] = 1.;

  for (size_t ii = 0; ii < cols; ++ii)
    permutations[ii] = ii;

  // compute (squared) column norms
  FieldVector<FieldType, cols> col_norms(0.);
  for (size_t rr = 0; rr < rows; ++rr)
    for (size_t cc = 0; cc < cols; ++cc)
      col_norms[cc] += std::pow(A[rr][cc], 2);

  FieldVector<FieldType, rows> w(0.);

  for (size_t jj = 0; jj < cols; ++jj) {

    // Pivoting
    // swap column jj and column with greatest norm
    auto max_it = std::max_element(col_norms.begin() + jj, col_norms.end());
    size_t max_index = std::distance(col_norms.begin(), max_it);
    if (max_index != jj) {
      for (size_t rr = 0; rr < rows; ++rr)
        std::swap(A[rr][jj], A[rr][max_index]);
      std::swap(col_norms[jj], col_norms[max_index]);
      std::swap(permutations[jj], permutations[max_index]);
    }

    // Matrix update
    // Reduction by householder matrix
    FieldType normx(0);
    for (size_t rr = jj; rr < rows; ++rr)
      normx += std::pow(A[rr][jj], 2);
    normx = std::sqrt(normx);

    if (XT::Common::FloatCmp::ne(normx, 0.)) {
      const auto s = -sign(A[jj][jj]);
      const FieldType u1 = A[jj][jj] - s * normx;
      w[jj] = 1.;
      for (size_t rr = jj + 1; rr < rows; ++rr)
        w[rr] = A[rr][jj] / u1;
      FieldType tau = -s * u1 / normx;
      // calculate A = H A
      multiply_householder_from_left(A, tau, w, jj, rows, jj, cols);
      multiply_householder_from_right(Q, tau, w, 0, rows, jj, cols);
    } // if (normx != 0)

    // Norm downdate
    for (size_t rr = jj + 1; rr < rows; ++rr)
      col_norms[rr] -= std::pow(A[jj][rr], 2);

  } // jj
} // void qr(...)


} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_ALGORITHMS_QR_HH
