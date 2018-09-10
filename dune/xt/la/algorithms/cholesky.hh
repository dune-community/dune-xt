// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017 - 2018)

#ifndef DUNE_XT_LA_ALGORITHMS_CHOLESKY_HH
#define DUNE_XT_LA_ALGORITHMS_CHOLESKY_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/math.hh>

#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/lapacke.hh>
#include <dune/xt/common/string.hh>

#include <dune/xt/la/container.hh>
#include <dune/xt/la/container/eye-matrix.hh>
#include <dune/xt/la/algorithms/triangular_solves.hh>

namespace Dune {
namespace XT {
namespace LA {
namespace internal {

// computes the LDL^T factorization of a tridiagonal matrix
template <class FirstVectorType, class SecondVectorType>
void tridiagonal_ldlt(FirstVectorType& diag, SecondVectorType& subdiag)
{
  using V1 = Common::VectorAbstraction<FirstVectorType>;
  using V2 = Common::VectorAbstraction<SecondVectorType>;
  V1::set_entry(diag, 0, V1::get_entry(diag, 0));
  for (size_t ii = 0; ii < diag.size() - 1; ++ii) {
    if (!(V1::get_entry(diag, ii) > 0)) // use !(.. > 0) instead of (.. <= 0) to also throw on NaNs
      DUNE_THROW(MathError, "LDL^T factorization failed!");
    V2::set_entry(subdiag, ii, V2::get_entry(subdiag, ii) / V1::get_entry(diag, ii));
    V1::set_entry(
        diag, ii + 1, V1::get_entry(diag, ii + 1) - V1::get_entry(diag, ii) * std::pow(V2::get_entry(subdiag, ii), 2));
  }
}

template <class FirstVectorType, class SecondVectorType, class VectorType>
std::enable_if_t<Common::is_vector<VectorType>::value, void>
solve_tridiag_ldlt(const FirstVectorType& diag, const SecondVectorType& subdiag, VectorType& vec)
{
  using V1 = Common::VectorAbstraction<FirstVectorType>;
  using V2 = Common::VectorAbstraction<SecondVectorType>;
  using V = Common::VectorAbstraction<VectorType>;
  using ScalarType = typename V::ScalarType;
  size_t size = vec.size();
  thread_local auto L =
      eye_matrix<CommonSparseMatrix<ScalarType>>(size, diagonal_pattern(size, size) + diagonal_pattern(size, size, -1));
  for (size_t ii = 0; ii < size - 1; ++ii)
    L.set_entry(ii + 1, ii, V2::get_entry(subdiag, ii));
  // solve LDL^T x = rhs;
  // first, solve L z = rhs;
  auto z = vec;
  solve_lower_triangular(L, z, vec);
  // solve D z = z
  for (size_t ii = 0; ii < size; ++ii)
    V::set_entry(z, ii, V::get_entry(z, ii) / V1::get_entry(diag, ii));
  // Now solve L^T x = z;
  solve_lower_triangular_transposed(L, vec, z);
}

template <class FirstVectorType, class SecondVectorType, class MatrixType>
std::enable_if_t<Common::is_matrix<MatrixType>::value, void>
solve_tridiag_ldlt(const FirstVectorType& diag, const SecondVectorType& subdiag, MatrixType& mat)
{
  using V1 = Common::VectorAbstraction<FirstVectorType>;
  using M = Common::MatrixAbstraction<MatrixType>;
  auto rhs_jj = diag;
  for (size_t jj = 0; jj < M::cols(mat); ++jj) {
    for (size_t ii = 0; ii < M::rows(mat); ++ii)
      V1::set_entry(rhs_jj, ii, M::get_entry(mat, ii, jj));
    solve_tridiag_ldlt(diag, subdiag, rhs_jj);
    for (size_t ii = 0; ii < M::rows(mat); ++ii)
      M::set_entry(mat, ii, jj, V1::get_entry(rhs_jj, ii));
  }
}

// computes the LL^T factorization of a symmetric positive definite matrix
template <bool only_set_nonzero, class MatrixType>
void cholesky_rowwise(MatrixType& A)
{
  using M = Common::MatrixAbstraction<MatrixType>;
  size_t size = M::rows(A);
  auto& L = A;
  for (size_t ii = 0; ii < size; ++ii) {
    for (size_t jj = 0; jj < ii; ++jj) {
      auto L_ij = M::get_entry(A, ii, jj);
      for (size_t kk = 0; kk < jj; ++kk)
        L_ij -= M::get_entry(L, ii, kk) * M::get_entry(L, jj, kk);
      L_ij /= M::get_entry(L, jj, jj);
      if (!only_set_nonzero || !XT::Common::is_zero(L_ij))
        M::set_entry(L, ii, jj, L_ij);
    } // jj
    auto L_ii = M::get_entry(A, ii, ii);
    for (size_t kk = 0; kk < ii; ++kk)
      L_ii -= std::pow(M::get_entry(L, ii, kk), 2);
    if (!(L_ii > 0)) // use !(.. > 0) instead of (.. <= 0) to also throw on NaNs
      DUNE_THROW(MathError, "Cholesky factorization failed!");
    M::set_entry(L, ii, ii, std::sqrt(L_ii));
  } // ii
}

template <bool only_set_nonzero, class MatrixType>
void cholesky_colwise(MatrixType& A)
{
  using M = Common::MatrixAbstraction<MatrixType>;
  size_t size = M::rows(A);
  auto& L = A;
  for (size_t jj = 0; jj < size; ++jj) {
    auto L_jj = M::get_entry(A, jj, jj);
    for (size_t kk = 0; kk < jj; ++kk)
      L_jj -= std::pow(M::get_entry(L, jj, kk), 2);
    if (!(L_jj > 0)) // use !(.. > 0) instead of (.. <= 0) to also throw on NaNs
      DUNE_THROW(MathError, "Cholesky factorization failed!");
    L_jj = std::sqrt(L_jj);
    M::set_entry(L, jj, jj, L_jj);
    const auto L_jj_inv = 1. / L_jj;
    for (size_t ii = jj + 1; ii < size; ++ii) {
      auto L_ij = M::get_entry(A, ii, jj);
      for (size_t kk = 0; kk < jj; ++kk)
        L_ij -= M::get_entry(L, ii, kk) * M::get_entry(L, jj, kk);
      L_ij *= L_jj_inv;
      if (!only_set_nonzero || !XT::Common::is_zero(L_ij))
        M::set_entry(L, ii, jj, L_ij);
    } // ii
  } // jj
}

template <class MatrixType>
typename std::enable_if_t<Common::MatrixAbstraction<MatrixType>::storage_layout == Common::StorageLayout::csr, void>
cholesky_csr(MatrixType& A)
{
  const auto* entries = A.entries();
  const auto* row_pointers = A.outer_index_ptr();
  const auto* column_indices = A.inner_index_ptr();
  using M = Common::MatrixAbstraction<MatrixType>;
  size_t size = M::rows(A);
  auto& L = A;
  for (size_t ii = 0; ii < size; ++ii) {
    for (size_t jj = 0; jj < ii; ++jj) {
      auto L_ij = M::get_entry(A, ii, jj);
      auto ll = row_pointers[ii];
      auto kk = row_pointers[jj];
      while (ll < row_pointers[ii + 1] && kk < row_pointers[jj + 1] && size_t(column_indices[kk]) < jj) {
        if (column_indices[ll] < column_indices[kk])
          ++ll;
        else if (column_indices[ll] > column_indices[kk])
          ++kk;
        else
          L_ij -= entries[ll++] * entries[kk++];
      }
      L_ij /= M::get_entry(L, jj, jj);
      if (!XT::Common::is_zero(L_ij))
        M::set_entry(L, ii, jj, L_ij);
    } // jj
    auto L_ii = M::get_entry(A, ii, ii);
    auto kk = row_pointers[ii];
    while (kk < row_pointers[ii + 1] && size_t(column_indices[kk]) < ii)
      L_ii -= std::pow(entries[kk++], 2);
    if (!(L_ii > 0)) // use !(.. > 0) instead of (.. <= 0) to also throw on NaNs
      DUNE_THROW(MathError, "Cholesky factorization failed!");
    M::set_entry(L, ii, ii, std::sqrt(L_ii));
  } // ii
}

template <class MatrixType>
typename std::enable_if_t<Common::MatrixAbstraction<MatrixType>::storage_layout == Common::StorageLayout::csc, void>
cholesky_csc(MatrixType& A)
{
  cholesky_colwise<true>(A);
}

template <class MatrixType>
typename std::enable_if_t<Common::MatrixAbstraction<MatrixType>::storage_layout != Common::StorageLayout::csr, void>
cholesky_csr(MatrixType& /*A*/)
{
  DUNE_THROW(Dune::InvalidStateException, "This only makes sense for matrices with compressed sparse row layout!");
}

template <class MatrixType>
typename std::enable_if_t<Common::MatrixAbstraction<MatrixType>::storage_layout != Common::StorageLayout::csc, void>
cholesky_csc(MatrixType& /*A*/)
{
  DUNE_THROW(Dune::InvalidStateException, "This only makes sense for matrices with compressed sparse column layout!");
}

template <class MatrixType,
          Common::StorageLayout storage_layout = Common::MatrixAbstraction<MatrixType>::storage_layout>
struct CholeskySolver
{
  using M = Common::MatrixAbstraction<MatrixType>;
  using ScalarType = typename M::ScalarType;
  static constexpr bool only_set_nonzero = storage_layout == Common::StorageLayout::csr
                                           || storage_layout == Common::StorageLayout::csc
                                           || std::is_base_of<IstlRowMajorSparseMatrix<ScalarType>, MatrixType>::value;

  static void cholesky(MatrixType& A)
  {
    const size_t size = M::rows(A);
    if (size != M::cols(A))
      DUNE_THROW(Dune::InvalidStateException, "Matrix has to be square!");
    if (false) {
#if HAVE_MKL || HAVE_LAPACKE
    } else if ((storage_layout == Common::StorageLayout::dense_row_major
                || storage_layout == Common::StorageLayout::dense_column_major)
               && size > 10) {
      const int lapacke_storage_layout = (storage_layout == Common::StorageLayout::dense_row_major)
                                             ? Common::Lapacke::row_major()
                                             : Common::Lapacke::col_major();
      assert(size <= std::numeric_limits<int>::max());
      int info = Common::Lapacke::dpotrf_work(
          lapacke_storage_layout, 'L', static_cast<int>(size), M::data(A), static_cast<int>(size));
      if (info)
        DUNE_THROW(Dune::MathError,
                   "Cholesky factorization using Lapacke::dpotrf failed with status " + Common::to_string(info) + "!");
#endif // HAVE_MKL || HAVE_LAPACKE
    } else if (storage_layout == Common::StorageLayout::csr)
      cholesky_csr(A);
    else if (storage_layout == Common::StorageLayout::csc)
      cholesky_csc(A);
    else
      cholesky_rowwise<only_set_nonzero>(A);
  } // static void cholesky(...)
}; // struct CholeskySolver<...>


template <class FirstVectorType, class SecondVectorType, class RhsType>
struct LDLTSolver
{
  using M = Common::MatrixAbstraction<RhsType>;
  using V = Common::VectorAbstraction<RhsType>;
  using V1 = Common::VectorAbstraction<FirstVectorType>;
  using V2 = Common::VectorAbstraction<SecondVectorType>;
  using ScalarType = typename std::conditional<M::is_matrix, typename M::ScalarType, typename V::ScalarType>::type;
  static constexpr bool is_row_major = M::is_matrix && M::storage_layout == Common::StorageLayout::dense_row_major;
  static constexpr bool is_col_major = (M::is_matrix && M::storage_layout == Common::StorageLayout::dense_column_major)
                                       || (V::is_vector && V::is_contiguous);
  static constexpr bool is_contiguous = is_row_major || is_col_major;
  static constexpr bool only_set_nonzero =
      M::is_matrix
      && (M::storage_layout == Common::StorageLayout::csr || M::storage_layout == Common::StorageLayout::csc
          || std::is_base_of<IstlRowMajorSparseMatrix<ScalarType>, RhsType>::value);

  static void tridiagonal_ldlt(FirstVectorType& diag, SecondVectorType& subdiag)
  {
    const size_t size = diag.size();
    if (subdiag.size() != size - 1)
      DUNE_THROW(InvalidStateException, "Wrong size of diag and subdiag!");
#if HAVE_MKL || HAVE_LAPACKE
    assert(size <= std::numeric_limits<int>::max());
    auto info = Common::Lapacke::dpttrf(static_cast<int>(size), Common::data(diag), Common::data(subdiag));
    if (info)
      DUNE_THROW(Dune::MathError, "Lapacke_dpptrf returned an error code!");
#else // HAVE_MKL || HAVE_LAPACKE
    internal::tridiagonal_ldlt(diag, subdiag);
#endif
  } // static void tridiagonal_ldlt(...)

  static void
  solve_tridiagonal_ldlt_factorized(const FirstVectorType& diag, const SecondVectorType& subdiag, RhsType& rhs)
  {
    const size_t size = diag.size();
    assert(subdiag.size() == size - 1);
    if (false) {
      ;
#if HAVE_MKL || HAVE_LAPACKE
    } else if (is_contiguous) {
      assert(V::is_vector || std::max(M::cols(rhs), size) <= std::numeric_limits<int>::max());
      int rhs_cols = V::is_vector ? 1 : int(M::cols(rhs));
      int info = Common::Lapacke::dpttrs(is_row_major ? Common::Lapacke::row_major() : Common::Lapacke::col_major(),
                                         static_cast<int>(size),
                                         rhs_cols,
                                         Common::data(diag),
                                         Common::data(subdiag),
                                         Common::data(rhs),
                                         is_row_major ? rhs_cols : static_cast<int>(size));
      if (info)
        DUNE_THROW(Dune::MathError, "Lapack dpttrs failed!");
#endif // HAVE_MKL || HAVE_LAPACKE
    } else {
      internal::solve_tridiag_ldlt(diag, subdiag, rhs);
    }

  } // static void solve_tridiagonal_ldlt_factorized(...)
}; // struct CholeskySolver<...>


} // namespace internal


template <class MatrixType>
typename std::enable_if_t<Common::is_matrix<MatrixType>::value, void> cholesky(MatrixType& A)
{
  internal::CholeskySolver<MatrixType>::cholesky(A);
} // void solve_lower_triangular(...)

template <class MatrixType, class VectorType>
typename std::enable_if_t<Common::is_matrix<MatrixType>::value && Common::is_vector<VectorType>::value, void>
solve_cholesky_factorized(const MatrixType& L, VectorType& rhs)
{
  auto x = rhs;
  solve_lower_triangular(L, x, rhs);
  solve_lower_triangular_transposed(L, rhs, x);
} // void solve_lower_triangular(...)


template <class FirstVectorType, class SecondVectorType>
typename std::enable_if_t<Common::is_vector<FirstVectorType>::value && Common::is_vector<SecondVectorType>::value, void>
tridiagonal_ldlt(FirstVectorType& diag, SecondVectorType& subdiag)
{
  internal::LDLTSolver<FirstVectorType, SecondVectorType, SecondVectorType>::tridiagonal_ldlt(diag, subdiag);
} // void tridiagonal_ldlt(...)


template <class FirstVectorType, class SecondVectorType, class RhsType>
typename std::enable_if_t<Common::is_vector<FirstVectorType>::value && Common::is_vector<SecondVectorType>::value
                              && (Common::is_vector<RhsType>::value || Common::is_matrix<RhsType>::value),
                          void>
solve_tridiagonal_ldlt_factorized(const FirstVectorType& diag, const SecondVectorType& subdiag, RhsType& rhs)
{
  internal::LDLTSolver<FirstVectorType, SecondVectorType, RhsType>::solve_tridiagonal_ldlt_factorized(
      diag, subdiag, rhs);
} // void tridiagonal_ldlt(...)


} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_ALGORITHMS_CHOLESKY_HH
