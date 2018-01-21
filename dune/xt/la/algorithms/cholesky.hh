// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_XT_LA_ALGORITHMS_CHOLESKY_HH
#define DUNE_XT_LA_ALGORITHMS_CHOLESKY_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/math.hh>

#include <dune/xt/common/float_cmp.hh>

#include <dune/xt/la/container/common/matrix/sparse.hh>
#include <dune/xt/la/container/common/vector/sparse.hh>

namespace Dune {
namespace XT {
namespace LA {


// copied and adapted from dune/geometry/affinegeometry.hh
template <class FirstMatrixImp, class SecondMatrixImp>
static bool cholesky_L(const Dune::DenseMatrix<FirstMatrixImp>& H, Dune::DenseMatrix<SecondMatrixImp>& L)
{
  typedef typename FirstMatrixImp::value_type FieldType;
  size_t size = H.size();
  for (int ii = 0; ii < size; ++ii) {
    FieldType& rii = L[ii][ii];

    FieldType xDiag = H[ii][ii];
    for (int jj = 0; jj < ii; ++jj)
      xDiag -= std::pow(L[ii][jj], 2);

    if (xDiag <= FieldType(0))
      return false;

    rii = std::sqrt(xDiag);

    FieldType invrii = FieldType(1) / rii;
    for (int ll = ii + 1; ll < size; ++ll) {
      FieldType x = H[ll][ii];
      for (int jj = 0; jj < ii; ++jj)
        x -= L[ii][jj] * L[ll][jj];
      L[ll][ii] = invrii * x;
    }
  }
  return true;
}

template <class DuneDenseMatrixImp, class MatrixTraits>
static bool cholesky_L(const Dune::DenseMatrix<DuneDenseMatrixImp>& H,
                       MatrixInterface<MatrixTraits, typename DuneDenseMatrixImp::value_type>& L)
{
  typedef typename DuneDenseMatrixImp::value_type FieldType;
  size_t size = H.size();
  for (int ii = 0; ii < size; ++ii) {
    FieldType xDiag = H[ii][ii];
    for (int jj = 0; jj < ii; ++jj)
      xDiag -= std::pow(L.get_entry(ii, jj), 2);

    if (xDiag <= FieldType(0))
      return false;

    L.set_entry(ii, ii, std::sqrt(xDiag));

    FieldType invrii = FieldType(1) / L.get_entry(ii, ii);
    for (int ll = ii + 1; ll < size; ++ll) {
      FieldType x = H[ll][ii];
      for (int jj = 0; jj < ii; ++jj)
        x -= L.get_entry(ii, jj) * L.get_entry(ll, jj);
      L.set_entry(ll, ii, invrii * x);
    }
  }
  return true;
}

template <class DuneDenseMatrixImp>
static bool cholesky_L(const Dune::DenseMatrix<DuneDenseMatrixImp>& H,
                       CommonSparseMatrixCsc<typename DuneDenseMatrixImp::value_type>& L)
{
  typedef typename DuneDenseMatrixImp::value_type FieldType;
  size_t size = H.size();
  thread_local std::vector<CommonSparseVector<FieldType>> rows;
  rows.resize(size, CommonSparseVector<FieldType>(size, size_t(0)));
  for (auto& vec : rows)
    vec.clear();
  L.clear();
  for (int ii = 0; ii < size; ++ii) {
    FieldType xDiag = H[ii][ii];
    for (size_t kk = 0; kk < rows[ii].entries().size(); ++kk)
      xDiag -= std::pow(rows[ii].entries()[kk], 2);

    if (XT::Common::FloatCmp::le(xDiag, FieldType(0)))
      return false;

    L.entries().push_back(std::sqrt(xDiag));
    L.row_indices().push_back(ii);
    rows[ii].set_new_entry(ii, L.entries().back());

    FieldType invrii = FieldType(1) / L.entries().back();
    for (int ll = ii + 1; ll < size; ++ll) {
      FieldType x = H[ll][ii] - rows[ii] * rows[ll];
      if (XT::Common::FloatCmp::ne(x, 0.)) {
        L.entries().push_back(invrii * x);
        L.row_indices().push_back(ll);
        rows[ll].set_new_entry(ii, L.entries().back());
      }
    } // ll
    L.column_pointers()[ii + 1] = L.row_indices().size();
  } // ii
  return true;
}

template <class DenseMatrixImp, class SparseMatrixImp, class DuneDenseMatrixImp>
static bool cholesky_L(const Dune::DenseMatrix<DuneDenseMatrixImp>& H,
                       CommonSparseOrDenseMatrix<DenseMatrixImp, SparseMatrixImp>& L)
{
  return L.sparse() ? cholesky_L(H, L.sparse_matrix()) : cholesky_L(H, L.dense_matrix());
}

template <class FieldType, int dim, int size>
static bool cholesky_L(const FieldVector<FieldMatrix<FieldType, dim, dim>, size>& H,
                       FieldVector<FieldMatrix<FieldType, dim, dim>, size>& L)
{
  for (size_t jj = 0; jj < size; ++jj) {
    bool chol_flag = cholesky_L(H[jj], L[jj]);
    if (!chol_flag)
      return false;
  }
  return true;
}

template <class FieldType, int size>
static bool cholesky_L(const FieldVector<FieldMatrix<FieldType, 2, 2>, size>& H,
                       FieldVector<FieldMatrix<FieldType, 2, 2>, size>& L)
{
  for (size_t jj = 0; jj < size; ++jj) {
    if (H[jj][0][0] <= 0.)
      return false;
    L[jj][0][0] = std::sqrt(H[jj][0][0]);
    L[jj][1][0] = H[jj][1][0] / L[jj][0][0];
    const auto val = H[jj][1][1] - std::pow(L[jj][1][0], 2);
    if (val <= 0.)
      return false;
    L[jj][1][1] = std::sqrt(val);
  }
  return true;
}

} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_ALGORITHMS_CHOLESKY_HH
