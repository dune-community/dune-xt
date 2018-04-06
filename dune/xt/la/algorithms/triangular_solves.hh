// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_XT_LA_ALGORITHMS_SOLVE_LOWER_TRIANGULAR_HH
#define DUNE_XT_LA_ALGORITHMS_SOLVE_LOWER_TRIANGULAR_HH

#include <dune/common/exceptions.hh>

#include <dune/xt/common/lapacke.hh>
#include <dune/xt/common/matrix.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/common/vector.hh>

#include <dune/xt/la/container/common/matrix/sparse.hh>

namespace Dune {
namespace XT {
namespace LA {
namespace internal {


// A simple forward solve for Ax = b with a lower triangular matrix A
// If transposed is true, A is assumed to be upper triangular and the
// equation A^T x = b is solved
template <class MatrixType, class VectorType, bool transposed>
void forward_solve(const MatrixType& A, VectorType& x, VectorType& rhs)
{
  typedef Common::MatrixAbstraction<MatrixType> M;
  for (size_t ii = 0; ii < M::rows(A); ++ii) {
    for (size_t jj = 0; jj < ii; ++jj)
      rhs[ii] -= M::get_entry(A, transposed ? jj : ii, transposed ? ii : jj) * x[jj];
    // We could simply use ==, but that gives a warning with -Wfloat-equal
    if (M::get_entry(A, ii, ii) == 0.)
      DUNE_THROW(Dune::MathError, "Triangular solve failed, matrix is singular!");
    x[ii] = rhs[ii] / M::get_entry(A, ii, ii);
  } // ii
}

// A simple backward solve for Ax = b with an upper triangular matrix A
// If transposed is true, A is assumed to be lower triangular and the
// equation A^T x = b is solved
template <class MatrixType, class VectorType, bool transposed>
void backward_solve(const MatrixType& A, VectorType& x, VectorType& rhs)
{
  typedef Common::MatrixAbstraction<MatrixType> M;
  const size_t num_rows = M::rows(A);
  const size_t num_cols = num_rows;
  for (int ii = static_cast<int>(num_rows) - 1; ii >= 0.; --ii) {
    for (size_t jj = ii + 1; jj < num_cols; ++jj)
      rhs[ii] -= M::get_entry(A, transposed ? jj : ii, transposed ? ii : jj) * x[jj];
    if (M::get_entry(A, ii, ii) == 0.)
      DUNE_THROW(Dune::MathError, "Triangular solve failed, matrix is singular!");
    x[ii] = rhs[ii] / M::get_entry(A, ii, ii);
  } // ii
}

template <class MatrixType>
constexpr bool has_compressed_sparse_layout = (Common::MatrixAbstraction<MatrixType>::storage_layout
                                               == Common::StorageLayout::csr)
                                              || (Common::MatrixAbstraction<MatrixType>::storage_layout
                                                  == Common::StorageLayout::csc);

// A forward solve to solve Ax = b with A lower triangular csr matrix
template <class MatrixType, class VectorType>
typename std::enable_if_t<has_compressed_sparse_layout<MatrixType>, void>
forward_solve_csr(const MatrixType& A, VectorType& x, VectorType& rhs)
{
  const auto* entries = A.entries();
  const auto* row_pointers = A.outer_index_ptr();
  const auto* column_indices = A.inner_index_ptr();
  size_t num_rows = A.rows();
  for (size_t ii = 0; ii < num_rows; ++ii) {
    // row_pointers[ii+1]-1 is the diagonal entry as we assume a lower triangular matrix with non-zero entries on the
    // diagonal
    auto kk = row_pointers[ii];
    for (; kk < row_pointers[ii + 1] - 1; ++kk)
      rhs[ii] -= entries[kk] * x[column_indices[kk]];
    if (entries[kk] == 0.)
      DUNE_THROW(Dune::MathError, "Triangular solve failed, matrix is singular!");
    x[ii] = rhs[ii] / entries[kk];
  } // ii
}

// A forward solve to solve Ax = b with A lower triangular csc matrix
template <class MatrixType, class VectorType>
typename std::enable_if_t<has_compressed_sparse_layout<MatrixType>, void>
forward_solve_csc(const MatrixType& A, VectorType& x, VectorType& rhs)
{
  const auto* entries = A.entries();
  const auto* column_pointers = A.outer_index_ptr();
  const auto* row_indices = A.inner_index_ptr();
  size_t num_rows = A.rows();
  for (size_t ii = 0; ii < num_rows; ++ii) {
    // column_pointers[ii] is the diagonal entry as we assume a lower triangular matrix with non-zero entries on the
    // diagonal
    auto kk = column_pointers[ii];
    if (entries[kk] == 0.)
      DUNE_THROW(Dune::MathError, "Triangular solve failed, matrix is singular!");
    x[ii] = rhs[ii] / entries[kk];
    ++kk;
    for (; kk < column_pointers[ii + 1]; ++kk)
      rhs[row_indices[kk]] -= entries[kk] * x[ii];
  } // ii
}

// A backward solve to solve Ax = b with A upper triangular csr matrix
template <class MatrixType, class VectorType>
typename std::enable_if_t<has_compressed_sparse_layout<MatrixType>, void>
backward_solve_csr(const MatrixType& A, VectorType& x, VectorType& rhs)
{
  const auto* entries = A.entries();
  const auto* row_pointers = A.outer_index_ptr();
  const auto* column_indices = A.inner_index_ptr();
  assert(A.rows() <= std::numeric_limits<int>::max());
  for (int ii = static_cast<int>(A.rows()) - 1; ii >= 0; ii--) {
    // row_pointers[ii] is the diagonal entry as we assume a upper triangular matrix with non-zero entries on the
    // diagonal
    const auto ll = row_pointers[ii];
    for (auto kk = ll + 1; kk < row_pointers[ii + 1]; ++kk)
      rhs[ii] -= entries[kk] * x[column_indices[kk]];
    if (entries[ll] == 0.)
      DUNE_THROW(Dune::MathError, "Triangular solve failed, matrix is singular!");
    x[ii] = rhs[ii] / entries[ll];
  } // ii
}

// A backward solve to solve Ax = b with A upper triangular csc matrix
template <class MatrixType, class VectorType>
typename std::enable_if_t<has_compressed_sparse_layout<MatrixType>, void>
backward_solve_csc(const MatrixType& A, VectorType& x, VectorType& rhs)
{
  const auto* entries = A.entries();
  const auto* column_pointers = A.outer_index_ptr();
  const auto* row_indices = A.inner_index_ptr();
  assert(A.rows() <= std::numeric_limits<int>::max());
  for (int ii = static_cast<int>(A.rows()) - 1; ii >= 0; ii--) {
    // column_pointers[ii+1]-1 is the diagonal entry as we assume an upper triangular matrix with non-zero entries
    // on the diagonal
    int kk = int(column_pointers[ii + 1]) - 1;
    if (entries[kk] == 0.)
      DUNE_THROW(Dune::MathError, "Triangular solve failed, matrix is singular!");
    x[ii] = rhs[ii] / entries[kk];
    --kk;
    for (; kk >= int(column_pointers[ii]); --kk)
      rhs[row_indices[kk]] -= entries[kk] * x[ii];
  } // ii
}

template <class MatrixType, class VectorType>
typename std::enable_if_t<!has_compressed_sparse_layout<MatrixType>, void>
backward_solve_csr(const MatrixType& /*A*/, VectorType& /*x*/, VectorType& /*rhs*/)
{
  DUNE_THROW(Dune::InvalidStateException,
             "This only makes sense for matrices with compressed sparse row/column layout!");
}

template <class MatrixType, class VectorType>
typename std::enable_if_t<!has_compressed_sparse_layout<MatrixType>, void>
forward_solve_csr(const MatrixType& /*A*/, VectorType& /*x*/, VectorType& /*rhs*/)
{
  DUNE_THROW(Dune::InvalidStateException,
             "This only makes sense for matrices with compressed sparse row/column layout!");
}

template <class MatrixType, class VectorType>
typename std::enable_if_t<!has_compressed_sparse_layout<MatrixType>, void>
backward_solve_csc(const MatrixType& /*A*/, VectorType& /*x*/, VectorType& /*rhs*/)
{
  DUNE_THROW(Dune::InvalidStateException,
             "This only makes sense for matrices with compressed sparse row/column layout!");
}

template <class MatrixType, class VectorType>
typename std::enable_if_t<!has_compressed_sparse_layout<MatrixType>, void>
forward_solve_csc(const MatrixType& /*A*/, VectorType& /*x*/, VectorType& /*rhs*/)
{
  DUNE_THROW(Dune::InvalidStateException,
             "This only makes sense for matrices with compressed sparse row/column layout!");
}


// TriangularSolver redirects to the appropiate forward or backward solve
template <class MatrixType,
          class VectorType,
          Common::MatrixPattern triangular_type,
          Common::Transpose transpose,
          Common::StorageLayout storage_layout = Common::MatrixAbstraction<MatrixType>::storage_layout>
struct TriangularSolver
{
  using M = typename Common::MatrixAbstraction<MatrixType>;
  using V = typename Common::VectorAbstraction<VectorType>;
  using ScalarType = typename M::ScalarType;

  static void solve(const MatrixType& A, VectorType& x)
  {
    auto* xp = V::data(x); // get pointer to x
    auto* rhs = xp; // use x to store rhs
    assert(std::max(M::rows(A), M::cols(A)) <= std::numeric_limits<int>::max());
    const int num_rows = static_cast<int>(M::rows(A));
    const int num_cols = static_cast<int>(M::cols(A));
    const bool trans = (transpose == Common::Transpose::yes);
    const bool lower = (triangular_type == Common::MatrixPattern::lower_triangular);
    const bool upper = !lower;
    const bool solve_forward = (lower && !trans) || (upper && trans);
    if (num_rows != num_cols)
      DUNE_THROW(Dune::InvalidStateException, "Matrix has to be square!");
    if (num_rows == 2) {
      if (solve_forward) {
        x[0] /= M::get_entry(A, 0, 0);
        x[1] = (x[1] - M::get_entry(A, trans ? 0 : 1, trans ? 1 : 0) * x[0]) / M::get_entry(A, 1, 1);
      } else {
        x[1] /= M::get_entry(A, 1, 1);
        x[0] = (x[0] - M::get_entry(A, trans ? 1 : 0, trans ? 0 : 1) * x[1]) / M::get_entry(A, 0, 0);
      } // if (solve_forward)
#if HAVE_MKL || HAVE_LAPACKE
    } else if (storage_layout == Common::StorageLayout::dense_row_major
               || storage_layout == Common::StorageLayout::dense_column_major) {
      const int blas_storage_layout = (storage_layout == Common::StorageLayout::dense_row_major)
                                          ? Common::Blas::row_major()
                                          : Common::Blas::col_major();
      const int blas_triangular =
          (triangular_type == Common::MatrixPattern::upper_triangular) ? Common::Blas::upper() : Common::Blas::lower();
      const int blas_trans = (transpose == Common::Transpose::yes) ? Common::Blas::trans() : Common::Blas::no_trans();
      trsm(blas_storage_layout,
           Common::Blas::left(),
           blas_triangular,
           blas_trans,
           Common::Blas::non_unit(),
           num_rows,
           1,
           ScalarType(1.),
           M::data(A),
           num_rows,
           rhs,
           storage_layout == Common::StorageLayout::dense_row_major ? 1 : num_rows);
#endif // HAVE_MKL || HAVE_LAPACKE
    } else if (storage_layout == Common::StorageLayout::csr) {
      if (lower && !trans) {
        forward_solve_csr(A, xp, rhs);
      } else if (lower && trans) {
        // A transposed lower triangular csr matrix is equivalent to an upper triangular csc matrix,
        // so use the backward_solve_csc
        backward_solve_csc(A, xp, rhs);
      } else if (upper && !trans) {
        backward_solve_csr(A, xp, rhs);
      } else { // (upper && trans)
        // A transposed upper triangular csr matrix is equivalent to a lower triangular csc matrix,
        // so use the forward_solve_csc
        forward_solve_csc(A, xp, rhs);
      }
    } else if (storage_layout == Common::StorageLayout::csc) {
      if (lower && !trans) {
        forward_solve_csc(A, xp, rhs);
      } else if (lower && trans) {
        // A transposed lower triangular csc matrix is equivalent to an upper triangular csr matrix,
        // so use the backward_solve_csr
        backward_solve_csr(A, xp, rhs);
      } else if (upper && !trans) {
        backward_solve_csc(A, xp, rhs);
      } else { // (upper && trans)
        // A transposed upper triangular csc matrix is equivalent to a lower triangular csr matrix,
        // so use the forward_solve_csr
        forward_solve_csr(A, xp, rhs);
      }
    } else {
      if (solve_forward)
        forward_solve<MatrixType, decltype(xp), trans>(A, xp, rhs);
      else
        backward_solve<MatrixType, decltype(xp), trans>(A, xp, rhs);
    }
  } // static void solve(...)

#if HAVE_LAPACKE || HAVE_MKL
private:
  template <class ScalarType>
  static std::enable_if_t<Common::is_arithmetic<ScalarType>::value, void> trsm(const int layout,
                                                                               const int side,
                                                                               const int uplo,
                                                                               const int transa,
                                                                               const int diag,
                                                                               const int m,
                                                                               const int n,
                                                                               const ScalarType alpha,
                                                                               const ScalarType* a,
                                                                               const int lda,
                                                                               ScalarType* b,
                                                                               const int ldb)
  {
    return Common::Blas::dtrsm(layout, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
  }

  template <class ScalarType>
  static std::enable_if_t<Common::is_complex<ScalarType>::value, void> trsm(const int layout,
                                                                            const int side,
                                                                            const int uplo,
                                                                            const int transa,
                                                                            const int diag,
                                                                            const int m,
                                                                            const int n,
                                                                            const ScalarType alpha,
                                                                            const ScalarType* a,
                                                                            const int lda,
                                                                            ScalarType* b,
                                                                            const int ldb)
  {
    return Common::Blas::ztrsm(layout,
                               side,
                               uplo,
                               transa,
                               diag,
                               m,
                               n,
                               static_cast<const void*>(&alpha),
                               static_cast<const void*>(a),
                               lda,
                               static_cast<void*>(b),
                               ldb);
  }
#endif // HAVE_LAPACKE || HAVE_MKL
};

// specialization for CommonSparseOrDenseMatrix
template <class DenseMatrixType,
          class SparseMatrixType,
          class VectorType,
          Common::MatrixPattern triangular_type,
          Common::Transpose transpose>
class TriangularSolver<CommonSparseOrDenseMatrix<DenseMatrixType, SparseMatrixType>,
                       VectorType,
                       triangular_type,
                       transpose,
                       Common::StorageLayout::other>
{
  typedef CommonSparseOrDenseMatrix<DenseMatrixType, SparseMatrixType> MatrixType;
  typedef Common::VectorAbstraction<MatrixType> M;
  typedef Common::VectorAbstraction<VectorType> V;

  static void solve(const MatrixType& A, VectorType& x)
  {
    A.sparse() ? TriangularSolver<SparseMatrixType, VectorType, triangular_type, transpose>::solve(A.sparse_matrix(), x)
               : TriangularSolver<DenseMatrixType, VectorType, triangular_type, transpose>::solve(A.dense_matrix(), x);
  }
}; // class TriangularSolver<CommonSparseOrDenseMatrix, ...>


// triangular_helper filters out sparse vectors.
// Triangular solves don't necessarily profit from sparse vectors, so we just take care that the resulting vector
// is stored in a sparse format so other algorithms can profit from the sparsity
template <class MatrixType,
          class FirstVectorType,
          class SecondVectorType,
          Common::MatrixPattern triangular_type,
          Common::Transpose transpose,
          bool sparse_vector = !Common::VectorAbstraction<FirstVectorType>::is_contiguous
                               || !Common::VectorAbstraction<SecondVectorType>::is_contiguous>
struct triangular_helper
{
  static_assert(triangular_type == Common::MatrixPattern::lower_triangular
                    || triangular_type == Common::MatrixPattern::upper_triangular,
                "The matrix has to be either upper or lower triangular!");

  static void solve(const MatrixType& A, FirstVectorType& x, const SecondVectorType& b)
  {
    std::copy_n(b.begin(), b.size(), x.begin()); // copy data to x
    TriangularSolver<MatrixType, FirstVectorType, triangular_type, transpose>::solve(A, x);
  }
};

template <class MatrixType,
          class FirstVectorType,
          class SecondVectorType,
          Common::MatrixPattern triangular_type,
          Common::Transpose transpose>
struct triangular_helper<MatrixType, FirstVectorType, SecondVectorType, triangular_type, transpose, true>
{
  static_assert(triangular_type == Common::MatrixPattern::lower_triangular
                    || triangular_type == Common::MatrixPattern::upper_triangular,
                "The matrix has to be either upper or lower triangular!");

  typedef Common::VectorAbstraction<FirstVectorType> V1;
  typedef Common::VectorAbstraction<SecondVectorType> V2;
  typedef typename V1::ScalarType ScalarType;

  static void solve(const MatrixType& A, FirstVectorType& x, const SecondVectorType& b)
  {
    // copy data from b to dense x vector
    const size_t size = V2::size(b);
    assert(V1::size(x) == size);
    std::vector<ScalarType> x_dense(size, 0.);
    for (size_t ii = 0; ii < size; ++ii)
      x_dense[ii] = V1::get_entry(b, ii);
    // call dense vector specialization
    TriangularSolver<MatrixType, std::vector<ScalarType>, triangular_type, transpose>::solve(A, x_dense);
    // set sparse entries
    x.clear();
    for (size_t ii = 0; ii < size; ++ii)
      if (Common::FloatCmp::ne(x_dense[ii], 0.))
        x.set_new_entry(ii, x_dense[ii]);
  }
};


} // namespace internal


/**
 * \brief solve A x = b, where A is lower triangular
 */
template <class MatrixType, class FirstVectorType, class SecondVectorType>
typename std::enable_if_t<Common::is_matrix<MatrixType>::value && Common::is_vector<FirstVectorType>::value
                              && Common::is_vector<SecondVectorType>::value,
                          void>
solve_lower_triangular(const MatrixType& A, FirstVectorType& x, const SecondVectorType& b)
{
  internal::triangular_helper<MatrixType,
                              FirstVectorType,
                              SecondVectorType,
                              Common::MatrixPattern::lower_triangular,
                              Common::Transpose::no>::solve(A, x, b);
} // void solve_lower_triangular(...)

/**
 * \brief solve A^T x = b, where A is lower triangular
 */
template <class MatrixType, class FirstVectorType, class SecondVectorType>
typename std::enable_if_t<Common::is_matrix<MatrixType>::value && Common::is_vector<FirstVectorType>::value
                              && Common::is_vector<SecondVectorType>::value,
                          void>
solve_lower_triangular_transposed(const MatrixType& A, FirstVectorType& x, const SecondVectorType& b)
{
  internal::triangular_helper<MatrixType,
                              FirstVectorType,
                              SecondVectorType,
                              Common::MatrixPattern::lower_triangular,
                              Common::Transpose::yes>::solve(A, x, b);
} // void solve_lower_triangular_transposed(...)

/**
 * \brief solve A x = b, where A is upper triangular
 */
template <class MatrixType, class FirstVectorType, class SecondVectorType>
typename std::enable_if_t<Common::is_matrix<MatrixType>::value && Common::is_vector<FirstVectorType>::value
                              && Common::is_vector<SecondVectorType>::value,
                          void>
solve_upper_triangular(const MatrixType& A, FirstVectorType& x, const SecondVectorType& b)
{
  internal::triangular_helper<MatrixType,
                              FirstVectorType,
                              SecondVectorType,
                              Common::MatrixPattern::upper_triangular,
                              Common::Transpose::no>::solve(A, x, b);
} // void solve_upper_triangular(...)

/**
 * \brief solve A^T x = b, where A is upper triangular
 */
template <class MatrixType, class FirstVectorType, class SecondVectorType>
typename std::enable_if_t<Common::is_matrix<MatrixType>::value && Common::is_vector<FirstVectorType>::value
                              && Common::is_vector<SecondVectorType>::value,
                          void>
solve_upper_triangular_transposed(const MatrixType& A, FirstVectorType& x, const SecondVectorType& b)
{
  internal::triangular_helper<MatrixType,
                              FirstVectorType,
                              SecondVectorType,
                              Common::MatrixPattern::upper_triangular,
                              Common::Transpose::yes>::solve(A, x, b);
} // void solve_upper_triangular_transposed(...)


} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_ALGORITHMS_SOLVE_LOWER_TRIANGULAR_HH
