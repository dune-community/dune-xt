// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_XT_LA_ALGORITHMS_TRIANGULAR_HH
#define DUNE_XT_LA_ALGORITHMS_TRIANGULAR_HH

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/xt/common/float_cmp.hh>

#include <dune/xt/la/container/common/matrix/sparse.hh>
#include <dune/xt/la/container/common/vector/sparse.hh>

namespace Dune {
namespace XT {
namespace LA {


/** Lower triangular solves
 * \brief solve A x = b, where A is lower triangular
 */

template <class MatrixImp, class VectorImp>
void solve_lower_triangular(const Dune::DenseMatrix<MatrixImp>& A,
                            Dune::DenseVector<VectorImp>& x,
                            const Dune::DenseVector<VectorImp>& b)
{
  const size_t num_rows = A.rows();
  // copy assignment operator does not work correctly for DenseVector,
  // so we need to cast it to the derived type first
  auto& rhs = static_cast<VectorImp&>(x); // use x to store rhs
  rhs = static_cast<const VectorImp&>(b); // copy data
  // forward solve
  for (size_t ii = 0; ii < num_rows; ++ii) {
    for (size_t jj = 0; jj < ii; ++jj)
      rhs[ii] -= A[ii][jj] * x[jj];
    x[ii] = rhs[ii] / A[ii][ii];
  }
} // void solve_lower_triangular(Dune::DenseMatrix, ...)

template <class MatrixImp, class ScalarType>
void solve_lower_triangular(const Dune::DenseMatrix<MatrixImp>& A,
                            CommonSparseVector<ScalarType>& x,
                            const CommonSparseVector<ScalarType>& b)
{
  x.clear();
  const size_t num_rows = A.rows();
  // copy rhs to dense vector to speed up random access
  std::vector<ScalarType> rhs(b.size(), 0.);
  for (size_t kk = 0; kk < b.entries().size(); ++kk)
    rhs[b.indices()[kk]] = b.entries()[kk];
  // forward solve
  for (size_t ii = 0; ii < num_rows; ++ii) {
    for (size_t jj = 0; jj < ii; ++jj)
      rhs[ii] -= A[ii][jj] * x[jj];
    const auto x_ii = rhs[ii] / A[ii][ii];
    if (XT::Common::FloatCmp::ne(x_ii, 0.))
      x.set_new_entry(ii, x_ii);
  } // ii
} // void solve_lower_triangular(Dune::DenseMatrix, ...)

template <class MatrixTraits, class ScalarType, class VectorImp>
void solve_lower_triangular(const MatrixInterface<MatrixTraits, ScalarType>& A,
                            Dune::DenseVector<VectorImp>& x,
                            const Dune::DenseVector<VectorImp>& b)
{
  const size_t num_rows = A.rows();
  // copy assignment operator does not work correctly for DenseVector,
  // so we need to cast it to the derived type first
  auto& rhs = static_cast<VectorImp&>(x); // use x to store rhs
  rhs = static_cast<const VectorImp&>(b); // copy data
  // forward solve
  for (size_t ii = 0; ii < num_rows; ++ii) {
    for (size_t jj = 0; jj < ii; ++jj)
      rhs[ii] -= A.get_entry(ii, jj) * x[jj];
    x[ii] = rhs[ii] / A.get_entry(ii, ii);
  }
} // void solve_lower_triangular(MatrixInterface<...>, ...)

template <class MatrixTraits, class ScalarType, class VectorImp>
void solve_lower_triangular(const MatrixInterface<MatrixTraits, ScalarType>& A,
                            CommonSparseVector<ScalarType>& x,
                            const Dune::DenseVector<VectorImp>& rhs)
{
  x.clear();
  const size_t num_rows = A.rows();
  // forward solve
  for (size_t ii = 0; ii < num_rows; ++ii) {
    ScalarType rhs_ii = rhs[ii];
    for (size_t jj = 0; jj < ii; ++jj)
      rhs_ii -= A.get_entry(ii, jj) * x[jj];
    x[ii] = rhs_ii / A.get_entry(ii, ii);
  }
} // void solve_lower_triangular(MatrixInterface<...>, ...)

template <class MatrixTraits, class ScalarType>
void solve_lower_triangular(const MatrixInterface<MatrixTraits, ScalarType>& A,
                            CommonSparseVector<ScalarType>& x,
                            const CommonSparseVector<ScalarType>& b)
{
  // copy rhs to dense vector
  Dune::DynamicVector<ScalarType> rhs(b.size(), 0.);
  for (size_t kk = 0; kk < b.entries().size(); ++kk)
    rhs[b.indices()[kk]] = b.entries()[kk];
  solve_lower_triangular(A, x, rhs);
} // void solve_lower_triangular(Dune::DenseMatrix, ...)

template <class ScalarType, class VectorImp>
void solve_lower_triangular(const CommonSparseMatrixCsr<ScalarType>& A,
                            Dune::DenseVector<VectorImp>& x,
                            const Dune::DenseVector<VectorImp>& b)
{
  // copy assignment operator does not work correctly for DenseVector,
  // so we need to cast it to the derived type first
  auto& rhs = static_cast<VectorImp&>(x); // use x to store rhs
  rhs = static_cast<const VectorImp&>(b); // copy data
  const auto& entries = A.entries();
  const auto& row_pointers = A.row_pointers();
  const auto& column_indices = A.column_indices();
  // forward solve
  for (size_t ii = 0; ii < A.rows(); ++ii) {
    // row_pointers[ii+1]-1 is the diagonal entry as we assume a lower triangular matrix with non-zero entries on the
    // diagonal
    size_t kk = row_pointers[ii];
    for (; kk < row_pointers[ii + 1] - 1; ++kk)
      rhs[ii] -= entries[kk] * x[column_indices[kk]];
    x[ii] = rhs[ii] / entries[kk];
  }
}

template <class ScalarType>
void solve_lower_triangular(const CommonSparseMatrixCsc<ScalarType>& A,
                            CommonSparseVector<ScalarType>& x,
                            const CommonSparseVector<ScalarType>& b_in)
{
  x.clear();
  std::vector<ScalarType> b(b_in.size(), 0.);
  for (size_t kk = 0; kk < b_in.entries().size(); ++kk)
    b[b_in.indices()[kk]] = b_in.entries()[kk];
  const auto& entries = A.entries();
  const auto& column_pointers = A.column_pointers();
  const auto& row_indices = A.row_indices();
  // forward solve
  const size_t size = b.size();
  for (size_t ii = 0; ii < size; ++ii) {
    // column_pointers[ii] is the diagonal entry as we assume a lower triangular matrix with non-zero entries on the
    // diagonal
    const auto rhs_ii = b[ii];
    if (Common::FloatCmp::ne(rhs_ii, 0.)) {
      size_t kk = column_pointers[ii];
      const auto x_ii = rhs_ii / entries[kk];
      x.set_new_entry(ii, x_ii);
      ++kk;
      for (; kk < column_pointers[ii + 1]; ++kk)
        b[row_indices[kk]] -= entries[kk] * x_ii;
    }
  } // ii
}

template <class ScalarType, int size>
void solve_lower_triangular(const CommonSparseMatrixCsc<ScalarType>& A,
                            CommonSparseVector<ScalarType>& x,
                            const FieldVector<ScalarType, size>& b_in)
{
  x.clear();
  auto b = b_in;
  const auto& entries = A.entries();
  const auto& column_pointers = A.column_pointers();
  const auto& row_indices = A.row_indices();
  // forward solve
  for (size_t ii = 0; ii < size; ++ii) {
    // column_pointers[ii] is the diagonal entry as we assume a lower triangular matrix with non-zero entries on the
    // diagonal
    const auto rhs_ii = b[ii];
    if (Common::FloatCmp::ne(rhs_ii, 0.)) {
      size_t kk = column_pointers[ii];
      const auto x_ii = rhs_ii / entries[kk];
      x.set_new_entry(ii, x_ii);
      ++kk;
      for (; kk < column_pointers[ii + 1]; ++kk)
        b[row_indices[kk]] -= entries[kk] * x_ii;
    }
  } // ii
}

template <class DenseMatrixImp, class SparseMatrixImp, class FirstVectorImp, class SecondVectorImp>
void solve_lower_triangular(const CommonSparseOrDenseMatrix<DenseMatrixImp, SparseMatrixImp>& A,
                            FirstVectorImp& x,
                            const SecondVectorImp& b)
{
  A.sparse() ? solve_lower_triangular(A.sparse_matrix(), x, b) : solve_lower_triangular(A.dense_matrix(), x, b);
} // void solve_lower_triangular(CommonSparseMatrixCsc, ...)

template <class FieldType>
void solve_lower_triangular(const Dune::FieldMatrix<FieldType, 2, 2>& A,
                            Dune::FieldVector<FieldType, 2>& x,
                            const Dune::FieldVector<FieldType, 2>& b)
{
  x[0] = b[0] / A[0][0];
  x[1] = (b[1] - A[1][0] * x[0]) / A[1][1];
} // void solve_lower_triangular


/** Upper triangular solves
 * \brief solve A x = b, where A is upper triangular
 */


template <class MatrixImp, class VectorImp>
void solve_upper_triangular(const Dune::DenseMatrix<MatrixImp>& A,
                            Dune::DenseVector<VectorImp>& x,
                            const Dune::DenseVector<VectorImp>& b)
{
  const size_t num_rows = A.rows();
  // copy assignment operator does not work correctly for DenseVector,
  // so we need to cast it to the derived type first
  auto& rhs = static_cast<VectorImp&>(x); // use x to store rhs
  rhs = static_cast<const VectorImp&>(b); // copy data
  const size_t num_cols = A.cols();
  // backward solve
  for (int ii = b.size() - 1; ii >= 0.; --ii) {
    for (size_t jj = ii + 1; jj < num_cols; ++jj)
      rhs[ii] -= A[ii][jj] * x[jj];
    x[ii] = rhs[ii] / A[ii][ii];
  }
} // void solve_upper_triangular(Dune::DenseMatrix, ...)

template <class MatrixTraits, class ScalarType, class VectorImp>
void solve_upper_triangular(const MatrixInterface<MatrixTraits, ScalarType>& A,
                            Dune::DenseVector<VectorImp>& x,
                            const Dune::DenseVector<VectorImp>& b)
{
  // copy assignment operator does not work correctly for DenseVector,
  // so we need to cast it to the derived type first
  auto& rhs = static_cast<VectorImp&>(x); // use x to store rhs
  rhs = static_cast<const VectorImp&>(b); // copy data
  const size_t num_cols = A.cols();
  // backward solve
  for (int ii = b.size() - 1; ii >= 0.; --ii) {
    for (size_t jj = ii + 1; jj < num_cols; ++jj)
      rhs[ii] -= A.get_entry(ii, jj) * x[jj];
    x[ii] = rhs[ii] / A.get_entry(ii, ii);
  }
} // void solve_upper_triangular(CommonDenseMatrix, ...)

template <class FieldType, int block_size, int num_intervals>
void solve_upper_triangular(const FieldVector<FieldMatrix<FieldType, block_size, block_size>, num_intervals>& A,
                            FieldVector<FieldType, block_size * num_intervals>& x,
                            const FieldVector<FieldType, block_size * num_intervals>& b)
{
  auto& rhs = x; // use x to store rhs
  rhs = b; // copy data
  // backward solve
  for (size_t jj = 0; jj < num_intervals; ++jj) {
    const auto offset = block_size * jj;
    for (int ll = block_size - 1; ll >= 0; --ll) {
      for (size_t mm = ll + 1; mm < block_size; ++mm)
        rhs[offset + ll] -= A[jj][ll][mm] * x[offset + mm];
      x[offset + ll] = rhs[offset + ll] / A[jj][ll][ll];
    } // ll
  } // jj
} // void solve_upper_triangular(CommonDenseMatrix, ...)

template <class ScalarType, class VectorImp>
void solve_upper_triangular(const CommonSparseMatrixCsc<ScalarType>& A,
                            Dune::DenseVector<VectorImp>& x,
                            const Dune::DenseVector<VectorImp>& b)
{
  assert(A.cols() == x.size());
  const auto& entries = A.entries();
  const auto& column_pointers = A.column_pointers();
  const auto& row_indices = A.row_indices();
  // copy assignment operator does not work correctly for DenseVector,
  // so we need to cast it to the derived type first
  auto& rhs = static_cast<VectorImp&>(x); // use x to store rhs
  rhs = static_cast<const VectorImp&>(b); // copy data
  // backsolve
  for (int ii = A.rows() - 1; ii >= 0; ii--) {
    // column_pointers[ii+1]-1 is the diagonal entry as we assume an upper triangular matrix with non-zero entries
    // on the diagonal
    int kk = int(column_pointers[ii + 1]) - 1;
    x[ii] = rhs[ii] / entries[kk];
    --kk;
    for (; kk >= int(column_pointers[ii]); --kk)
      rhs[row_indices[kk]] -= entries[kk] * x[ii];
  } // ii
} // void solve_upper_triangular(CommonSparseMatrixCsc, ...)

template <class DenseMatrixImp, class SparseMatrixImp, class FirstVectorImp, class SecondVectorImp>
void solve_upper_triangular(const CommonSparseOrDenseMatrix<DenseMatrixImp, SparseMatrixImp>& A,
                            FirstVectorImp& x,
                            const SecondVectorImp& b)
{
  A.sparse() ? solve_upper_triangular(A.sparse_matrix(), x, b) : solve_upper_triangular(A.dense_matrix(), x, b);
} // void solve_upper_triangular(CommonSparseMatrixCsc, ...)

/** Lower triangular transposed solves
 * \brief solve A^T x = b, where A is lower triangular
 */

// solve A^T x = b, where A is lower triangular
template <class MatrixImp, class FirstVectorImp, class SecondVectorImp>
void solve_lower_triangular_transposed(const Dune::DenseMatrix<MatrixImp>& A,
                                       Dune::DenseVector<FirstVectorImp>& x,
                                       const Dune::DenseVector<SecondVectorImp>& b)
{
  // copy assignment operator does not work correctly for DenseVector,
  // so we need to cast it to the derived type first
  auto& rhs = static_cast<FirstVectorImp&>(x); // use x to store rhs
  rhs = static_cast<const SecondVectorImp&>(b); // copy data
  // backsolve
  for (int ii = int(A.N()) - 1; ii >= 0; ii--) {
    for (size_t jj = ii + 1; jj < A.N(); jj++)
      rhs[ii] -= A[jj][ii] * x[jj];
    x[ii] = rhs[ii] / A[ii][ii];
  }
}

template <class MatrixTraits, class FirstVectorImp, class SecondVectorImp>
void solve_lower_triangular_transposed(const MatrixInterface<MatrixTraits, typename MatrixTraits::ScalarType>& A,
                                       Dune::DenseVector<FirstVectorImp>& x,
                                       const Dune::DenseVector<SecondVectorImp>& b)
{
  // copy assignment operator does not work correctly for DenseVector,
  // so we need to cast it to the derived type first
  auto& rhs = static_cast<FirstVectorImp&>(x); // use x to store rhs
  rhs = static_cast<const SecondVectorImp&>(b); // copy data
  // backsolve
  for (int ii = int(A.cols()) - 1; ii >= 0; ii--) {
    for (size_t jj = ii + 1; jj < A.rows(); jj++)
      rhs[ii] -= A.get_entry(jj, ii) * x[jj];
    x[ii] = rhs[ii] / A.get_entry(ii, ii);
  }
}

template <class MatrixTraits, class FirstVectorImp, class SecondVectorTraits>
void solve_lower_triangular_transposed(const MatrixInterface<MatrixTraits, typename MatrixTraits::ScalarType>& A,
                                       Dune::DenseVector<FirstVectorImp>& x,
                                       const VectorInterface<SecondVectorTraits, typename MatrixTraits::ScalarType>& b)
{
  // backsolve
  for (int ii = int(A.cols()) - 1; ii >= 0; ii--) {
    auto rhs_ii = b.get_entry(ii);
    for (size_t jj = ii + 1; jj < A.rows(); jj++)
      rhs_ii -= A.get_entry(jj, ii) * x[jj];
    x[ii] = rhs_ii / A.get_entry(ii, ii);
  }
}

template <class MatrixTraits, class VectorTraits>
void solve_lower_triangular_transposed(const MatrixInterface<MatrixTraits, typename MatrixTraits::ScalarType>& A,
                                       CommonSparseVector<typename MatrixTraits::ScalarType>& x,
                                       const VectorInterface<VectorTraits, typename MatrixTraits::ScalarType>& b)
{
  x.clear();
  // backsolve
  typename MatrixTraits::ScalarType rhs_ii;
  for (int ii = int(A.cols()) - 1; ii >= 0; --ii) {
    rhs_ii = b.get_entry(ii);
    for (size_t jj = ii + 1; jj < A.rows(); jj++)
      rhs_ii -= A.get_entry(jj, ii) * x.get_entry(jj);
    if (Common::FloatCmp::ne(rhs_ii, 0.))
      x.set_new_entry(ii, rhs_ii / A.get_entry(ii, ii), true);
  } // ii
}

template <class ScalarType, class FirstVectorImp, class SecondVectorImp>
void solve_lower_triangular_transposed(const CommonSparseMatrixCsr<ScalarType>& A,
                                       Dune::DenseVector<FirstVectorImp>& x,
                                       const Dune::DenseVector<SecondVectorImp>& b)
{
  // copy assignment operator does not work correctly for DenseVector,
  // so we need to cast it to the derived type first
  auto& rhs = static_cast<FirstVectorImp&>(x); // use x to store rhs
  rhs = static_cast<const SecondVectorImp&>(b); // copy data
  // backsolve
  for (int ii = int(A.rows()) - 1; ii >= 0; ii--) {
    for (size_t jj = ii + 1; jj < A.rows(); jj++)
      rhs[ii] -= A.get_entry(jj, ii) * x[jj];
    x[ii] = rhs[ii] / A.get_entry(ii, ii);
  }
}

template <class ScalarType, class VectorTraitsImp>
void solve_lower_triangular_transposed(const CommonSparseMatrixCsc<ScalarType>& A,
                                       CommonSparseVector<ScalarType>& x,
                                       const VectorInterface<VectorTraitsImp, ScalarType>& b)
{
  x.clear();
  const auto& entries = A.entries();
  const auto& column_pointers = A.column_pointers();
  const auto& row_indices = A.row_indices();
  // backsolve
  ScalarType rhs_ii;
  for (int ii = int(A.cols()) - 1; ii >= 0; --ii) {
    // column_pointers[ii] is the diagonal entry as we assume a lower triangular matrix with non-zero entries on the
    // diagonal
    rhs_ii = b.get_entry(ii);
    const size_t ll = column_pointers[ii];
    for (size_t kk = ll + 1; kk < column_pointers[ii + 1]; ++kk)
      rhs_ii -= entries[kk] * x.get_entry(row_indices[kk]);
    if (Common::FloatCmp::ne(rhs_ii, 0.))
      x.set_new_entry(ii, rhs_ii / entries[ll], true);
  } // ii
}

template <class ScalarType, int size>
void solve_lower_triangular_transposed(const CommonSparseMatrixCsc<ScalarType>& A,
                                       FieldVector<ScalarType, size>& x,
                                       const CommonSparseVector<ScalarType>& b)
{
  std::fill(x.begin(), x.end(), 0.);
  const auto& entries = A.entries();
  const auto& column_pointers = A.column_pointers();
  const auto& row_indices = A.row_indices();
  // backsolve
  ScalarType rhs_ii;
  for (int ii = int(A.cols()) - 1; ii >= 0; --ii) {
    // column_pointers[ii] is the diagonal entry as we assume a lower triangular matrix with non-zero entries on the
    // diagonal
    rhs_ii = b.get_entry(ii);
    const size_t ll = column_pointers[ii];
    for (size_t kk = ll + 1; kk < column_pointers[ii + 1]; ++kk)
      rhs_ii -= entries[kk] * x[row_indices[kk]];
    x[ii] = rhs_ii / entries[ll];
  } // ii
}

template <class DenseMatrixImp, class SparseMatrixImp, class FirstVectorImp, class SecondVectorImp>
void solve_lower_triangular_transposed(const CommonSparseOrDenseMatrix<DenseMatrixImp, SparseMatrixImp>& A,
                                       FirstVectorImp& x,
                                       const SecondVectorImp& b)
{
  A.sparse() ? solve_lower_triangular_transposed(A.sparse_matrix(), x, b)
             : solve_lower_triangular_transposed(A.dense_matrix(), x, b);
} // void solve_lower_triangular_transposed(CommonSparseMatrixCsc, ...)

template <class FieldType>
void solve_lower_triangular_transposed(const Dune::FieldMatrix<FieldType, 2, 2>& A,
                                       Dune::FieldVector<FieldType, 2>& x,
                                       const Dune::FieldVector<FieldType, 2>& b)
{
  x[1] = b[1] / A[1][1];
  x[0] = (b[0] - A[1][0] * x[1]) / A[0][0];
} // void solve_lower_triangular


} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_ALGORITHMS_TRIANGULAR_HH
