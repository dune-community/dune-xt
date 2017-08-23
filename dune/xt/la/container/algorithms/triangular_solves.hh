// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2017 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2017)

#ifndef DUNE_XT_LA_CONTAINER_ALGORITHMS_TRIANGULAR_HH
#define DUNE_XT_LA_CONTAINER_ALGORITHMS_TRIANGULAR_HH

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

template <class ScalarType, int size>
void solve_lower_triangular(const FieldMatrix<ScalarType, size, size>& A,
                            FieldVector<ScalarType, size>& x,
                            const FieldVector<ScalarType, size>& b)
{
  FieldVector<ScalarType, size>& rhs = x; // use x to store rhs
  rhs = b; // copy data
  // forward solve
  for (size_t ii = 0; ii < size; ++ii) {
    for (size_t jj = 0; jj < ii; ++jj)
      rhs[ii] -= A[ii][jj] * x[jj];
    x[ii] = rhs[ii] / A[ii][ii];
  }
} // void solve_lower_triangular(FieldMatrix, ...)

template <class ScalarType, int size>
void solve_lower_triangular(const CommonSparseMatrixCsr<ScalarType>& A,
                            FieldVector<ScalarType, size>& x,
                            const FieldVector<ScalarType, size>& b)
{
  FieldVector<ScalarType, size>& rhs = x; // use x to store rhs
  rhs = b; // copy data
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
  static thread_local std::vector<ScalarType> b(b_in.size());
  std::fill(b.begin(), b.end(), 0.);
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

/** Upper triangular solves
 * \brief solve A x = b, where A is upper triangular
 */

template <class ScalarType, int size>
void solve_upper_triangular(const CommonSparseMatrixCsc<ScalarType>& A,
                            FieldVector<ScalarType, size>& x,
                            const FieldVector<ScalarType, size>& b)
{
  assert(A.cols() == size);
  const auto& entries = A.entries();
  const auto& column_pointers = A.column_pointers();
  const auto& row_indices = A.row_indices();
  FieldVector<ScalarType, size>& rhs = x; // use x to store rhs
  rhs = b; // copy data
  // backsolve
  for (int ii = size - 1; ii >= 0; ii--) {
    // column_pointers[ii+1]-1 is the diagonal entry as we assume an upper triangular matrix with non-zero entries
    // on the diagonal
    if (Common::FloatCmp::ne(rhs[ii], 0.)) {
      int kk = int(column_pointers[ii + 1]) - 1;
      x[ii] = rhs[ii] / entries[kk];
      --kk;
      for (; kk >= int(column_pointers[ii]); --kk)
        rhs[row_indices[kk]] -= entries[kk] * x[ii];
    }
  } // ii
} // void solve_upper_triangular(CommonSparseMatrixCsc, ...)


/** Lower triangular transposed solves
 * \brief solve A^T x = b, where A is lower triangular
 */

// solve A^T x = b, where A is lower triangular
template <class ScalarType, int size>
void solve_lower_triangular_transposed(const FieldMatrix<ScalarType, size, size>& A,
                                       FieldVector<ScalarType, size>& x,
                                       const FieldVector<ScalarType, size>& b)
{
  FieldVector<ScalarType, size>& rhs = x; // use x to store rhs
  rhs = b; // copy data
  // backsolve
  double min_eigval(std::abs(A[0][0]));
  double max_eigval = min_eigval;
  for (int ii = int(A.N()) - 1; ii >= 0; ii--) {
    auto abs = std::abs(A[ii][ii]);
    min_eigval = std::min(abs, min_eigval);
    max_eigval = std::max(abs, max_eigval);
    for (size_t jj = ii + 1; jj < A.N(); jj++)
      rhs[ii] -= A[jj][ii] * x[jj];
    x[ii] = rhs[ii] / A[ii][ii];
  }
  if (Common::FloatCmp::eq(min_eigval, 0.) || max_eigval / min_eigval > 1e10)
    DUNE_THROW(Dune::FMatrixError, "A is singular!");
}

template <class ScalarType, int size>
void solve_lower_triangular_transposed(const CommonSparseMatrixCsr<ScalarType>& A,
                                       FieldVector<ScalarType, size>& x,
                                       const FieldVector<ScalarType, size>& b)
{
  FieldVector<ScalarType, size>& rhs = x; // use x to store rhs
  rhs = b; // copy data
  // backsolve
  double min_eigval(std::abs(A.get_entry(0, 0)));
  double max_eigval = min_eigval;
  for (int ii = int(A.rows()) - 1; ii >= 0; ii--) {
    auto abs = std::abs(A.get_entry(ii, ii));
    min_eigval = std::min(abs, min_eigval);
    max_eigval = std::max(abs, max_eigval);
    for (size_t jj = ii + 1; jj < A.rows(); jj++)
      rhs[ii] -= A.get_entry(jj, ii) * x[jj];
    x[ii] = rhs[ii] / A.get_entry(ii, ii);
  }
  if (Common::FloatCmp::eq(min_eigval, 0.) || max_eigval / min_eigval > 1e10)
    DUNE_THROW(Dune::FMatrixError, "A is singular!");
}

template <class ScalarType>
void solve_lower_triangular_transposed(const CommonSparseMatrixCsc<ScalarType>& A,
                                       CommonSparseVector<ScalarType>& x,
                                       const CommonSparseVector<ScalarType>& b)
{
  x.clear();
  const auto& entries = A.entries();
  const auto& column_pointers = A.column_pointers();
  const auto& row_indices = A.row_indices();
  // forward solve
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
  // forward solve
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


} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_CONTAINER_ALGORITHMS_TRIANGULAR_HH
