// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2017 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Barbara Verfürth (2015)
//   Felix Schindler  (2013 - 2017)
//   Rene Milk        (2014 - 2016)
//   Tobias Leibner   (2014, 2016 - 2017)

#ifndef DUNE_XT_LA_CONTAINER_ALGORITHMS_CHOLESKY_HH
#define DUNE_XT_LA_CONTAINER_ALGORITHMS_CHOLESKY_HH

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
template <class FieldType, int size>
static bool cholesky_L(const FieldMatrix<FieldType, size, size>& H, FieldMatrix<FieldType, size, size>& L)
{
  for (int ii = 0; ii < size; ++ii) {
    FieldType& rii = L[ii][ii];

    FieldType xDiag = H[ii][ii];
    for (int jj = 0; jj < ii; ++jj)
      xDiag -= std::pow(L[ii][jj], 2);

    if (XT::Common::FloatCmp::le(xDiag, FieldType(0)))
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

template <class FieldType, int size>
static bool cholesky_L(const FieldMatrix<FieldType, size, size>& H, CommonSparseMatrixCsc<FieldType>& L)
{
  thread_local FieldVector<CommonSparseVector<FieldType>, size> rows(CommonSparseVector<FieldType>(size, size_t(0)));
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


} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_CONTAINER_ALGORITHMS_CHOLESKY_HH
