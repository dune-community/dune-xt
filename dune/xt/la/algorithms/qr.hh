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

#include <dune/xt/common/lapacke.hh>
#include <dune/xt/common/matrix.hh>
#include <dune/xt/common/vector.hh>

#include <dune/xt/la/algorithms/triangular_solves.hh>
#include <dune/xt/la/container/eye-matrix.hh>

namespace Dune {
namespace XT {
namespace LA {
namespace internal {


template <class FirstVectorType, class SecondVectorType>
typename Common::VectorAbstraction<FirstVectorType>::ScalarType
partial_dot(const FirstVectorType& a, const SecondVectorType& b, const size_t begin, const size_t end)
{
  typedef Common::VectorAbstraction<FirstVectorType> V1;
  typedef Common::VectorAbstraction<SecondVectorType> V2;
  typename V1::ScalarType ret(0);
  for (size_t ii = begin; ii < end; ++ii)
    ret += V1::get_entry(a, ii) * V2::get_entry(b, ii);
  return ret;
}

// Calculates x(begin:end) = H * x(begin:end) where H = I-tau*w*w^T and w = v[row_begin:row_end]
template <class FirstVectorType, class VectorType>
void multiply_householder_from_left(FirstVectorType& x,
                                    const typename Common::VectorAbstraction<VectorType>::ScalarType& tau,
                                    const VectorType& v,
                                    const size_t begin,
                                    const size_t end)
{
  typedef Common::VectorAbstraction<FirstVectorType> V1;
  typedef Common::VectorAbstraction<VectorType> V2;
  // calculate w^T x first
  auto wT_x = partial_dot(x, v, begin, end);
  for (size_t rr = begin; rr < end; ++rr)
    V1::add_to_entry(x, rr, -tau * wT_x * V2::get_entry(v, rr));
}

// Calculates A(row_begin:row_end,col_begin:col_end) = H * A(row_begin:row_end,col_begin:col_end) where H = I-tau*w*w^T
// and w = v[row_begin:row_end]
template <class MatrixType, class VectorType>
void multiply_householder_from_left(MatrixType& A,
                                    const typename Common::MatrixAbstraction<MatrixType>::ScalarType& tau,
                                    const VectorType& v,
                                    const size_t row_begin,
                                    const size_t row_end,
                                    const size_t col_begin,
                                    const size_t col_end)
{
  typedef Common::MatrixAbstraction<MatrixType> M;
  typedef Common::VectorAbstraction<VectorType> V;
  // calculate w^T A first
  VectorType wT_A = V::create(M::cols(A), 0.);
  for (size_t rr = row_begin; rr < row_end; ++rr)
    for (size_t cc = col_begin; cc < col_end; ++cc)
      V::add_to_entry(wT_A, cc, V::get_entry(v, rr) * M::get_entry(A, rr, cc));
  for (size_t rr = row_begin; rr < row_end; ++rr)
    for (size_t cc = col_begin; cc < col_end; ++cc)
      M::add_to_entry(A, rr, cc, -tau * V::get_entry(v, rr) * V::get_entry(wT_A, cc));
}

// Calculates A = A * H.
// \see multiply_householder_from_left
template <class MatrixType, class VectorType>
void multiply_householder_from_right(MatrixType& A,
                                     const typename Common::MatrixAbstraction<MatrixType>::ScalarType& tau,
                                     const VectorType& v,
                                     const size_t row_begin,
                                     const size_t row_end,
                                     const size_t col_begin,
                                     const size_t col_end)
{
  typedef Common::MatrixAbstraction<MatrixType> M;
  typedef Common::VectorAbstraction<VectorType> V;
  // calculate A w first
  VectorType Aw = V::create(M::rows(A), 0.);
  for (size_t rr = row_begin; rr < row_end; ++rr)
    for (size_t cc = col_begin; cc < col_end; ++cc)
      V::add_to_entry(Aw, cc, V::get_entry(v, cc) * M::get_entry(A, rr, cc));
  for (size_t rr = row_begin; rr < row_end; ++rr)
    for (size_t cc = col_begin; cc < col_end; ++cc)
      M::add_to_entry(A, rr, cc, -tau * V::get_entry(v, cc) * V::get_entry(Aw, rr));
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
template <class MatrixType, class VectorType, class IndexVectorType>
void qr_decomposition(MatrixType& A, VectorType& tau, IndexVectorType& permutations)
{
  typedef Common::MatrixAbstraction<MatrixType> M;
  typedef Common::VectorAbstraction<VectorType> V;
  typedef Common::VectorAbstraction<IndexVectorType> VI;
  typedef typename M::RealType RealType;

  const size_t num_rows = M::rows(A);
  const size_t num_cols = M::cols(A);
  assert(tau.size() == num_cols && permutations.size() == num_cols);
  std::fill(tau.begin(), tau.end(), 0.);
  for (size_t ii = 0; ii < num_cols; ++ii)
    VI::set_entry(permutations, ii, ii);

  // compute (squared) column norms
  std::vector<RealType> col_norms(num_cols);
  for (size_t rr = 0; rr < num_rows; ++rr)
    for (size_t cc = 0; cc < num_cols; ++cc)
      col_norms[cc] += std::pow(M::get_entry(A, rr, cc), 2);

  VectorType w = V::create(num_rows, 0.);

  for (size_t jj = 0; jj < num_cols - 1; ++jj) {

    // Pivoting
    // swap column jj and column with greatest norm
    auto max_it = std::max_element(col_norms.begin() + jj, col_norms.end());
    size_t max_index = std::distance(col_norms.begin(), max_it);
    if (max_index != jj) {
      auto tmp = col_norms[jj];
      col_norms[jj] = col_norms[max_index];
      col_norms[max_index] = tmp;
      auto tmp_index = VI::get_entry(permutations, jj);
      VI::set_entry(permutations, jj, VI::get_entry(permutations, max_index));
      VI::set_entry(permutations, max_index, tmp_index);
      for (size_t rr = 0; rr < num_rows; ++rr) {
        tmp = M::get_entry(A, rr, max_index);
        M::set_entry(A, rr, max_index, M::get_entry(A, rr, jj));
        M::set_entry(A, rr, jj, tmp);
      }
    }

    // Matrix update
    // Reduction by householder matrix
    auto normx = 0.;
    for (size_t rr = jj; rr < num_rows; ++rr)
      normx += std::pow(M::get_entry(A, rr, jj), 2);
    normx = std::sqrt(normx);

    if (normx != 0.) {
      const auto s = -sign(M::get_entry(A, jj, jj));
      const auto u1 = M::get_entry(A, jj, jj) - s * normx;
      V::set_entry(w, jj, 1.);
      for (size_t rr = jj + 1; rr < num_rows; ++rr) {
        V::set_entry(w, rr, M::get_entry(A, rr, jj) / u1);
        M::set_entry(A, rr, jj, V::get_entry(w, rr));
      }
      M::set_entry(A, jj, jj, s * normx);
      V::set_entry(tau, jj, -s * u1 / normx);
      // calculate A = H A
      multiply_householder_from_left(A, V::get_entry(tau, jj), w, jj, num_rows, jj + 1, num_cols);
    } // if (normx != 0)

    // Norm downdate
    for (size_t cc = jj + 1; cc < num_cols; ++cc)
      col_norms[cc] -= std::pow(M::get_entry(A, jj, cc), 2);

  } // jj
} // void qr_decomposition(...)

template <class MatrixType,
          class VectorType,
          class IndexVectorType = std::vector<int>,
          Common::StorageLayout storage_layout = Common::MatrixAbstraction<MatrixType>::storage_layout>
struct QrHelper
{
  typedef Common::MatrixAbstraction<MatrixType> M;
  typedef Common::VectorAbstraction<VectorType> V;
  typedef Common::VectorAbstraction<IndexVectorType> VI;
  static const bool is_row_major = (storage_layout == Common::StorageLayout::dense_row_major);
  static const bool has_contiguous_storage = (storage_layout == Common::StorageLayout::dense_row_major)
                                             || (storage_layout == Common::StorageLayout::dense_column_major);

  static int lapacke_storage_layout()
  {
    return storage_layout == Common::StorageLayout::dense_row_major ? Common::Lapacke::row_major()
                                                                    : Common::Lapacke::col_major();
  }

  static void qr(MatrixType& A, VectorType& tau, IndexVectorType& permutations)
  {
    if (false) {
      ;
#if HAVE_LAPACKE || HAVE_MKL
    } else if (has_contiguous_storage) {
      auto info = Common::Lapacke::dgeqp3(lapacke_storage_layout(),
                                          M::rows(A),
                                          M::cols(A),
                                          M::data(A),
                                          is_row_major ? M::cols(A) : M::rows(A),
                                          VI::data(permutations),
                                          V::data(tau));
      if (info)
        DUNE_THROW(Dune::MathError, "QR factorization failed");
      // Lapack indices are 1-based, convert to 0-based
      for (auto& index : permutations)
        index -= 1;
#endif // HAVE_LAPACKE || HAVE_MKL
    } else {
      qr_decomposition(A, tau, permutations);
    }
  } // static void qr(...)

  static typename M::template MatrixTypeTemplate<M::static_rows, M::static_rows>
  calculate_q_from_qr(const MatrixType& QR, const VectorType& tau)
  {
    if (M::cols(QR) > M::rows(QR))
      DUNE_THROW(NotImplemented, "Not yet implemented for cols > rows");
    if (false) {
      ;
#if HAVE_LAPACKE || HAVE_MKL
    } else if (has_contiguous_storage) {
      // create Q and copy values of QR to Q;
      const size_t num_rows = M::rows(QR);
      const size_t num_cols = M::cols(QR);
      auto ret = M::template create<M::static_rows, M::static_rows>(num_rows, num_rows, 0.);
      typedef Common::MatrixAbstraction<decltype(ret)> Mret;
      for (size_t ii = 0; ii < num_rows; ++ii)
        for (size_t jj = 0; jj < num_cols; ++jj)
          Mret::set_entry(ret, ii, jj, M::get_entry(QR, ii, jj));
      auto info = Common::Lapacke::dorgqr(
          lapacke_storage_layout(), num_rows, num_rows, num_cols, Mret::data(ret), num_rows, V::data(tau));
      if (info)
        DUNE_THROW(Dune::MathError, "Calculating Q explicitly failed!");
      return ret;
#endif // HAVE_LAPACKE || HAVE_MKL
    } else {
      return calculate_q_explicitly(QR, tau);
    }
  } // static void qr(...)

  template <XT::Common::Transpose transpose = XT::Common::Transpose::no,
            class SecondVectorType = VectorType,
            class ThirdVectorType = VectorType>
  static void
  apply_q_from_qr(const MatrixType& QR, const VectorType& tau, const SecondVectorType& x, ThirdVectorType& y)
  {
    typedef Common::VectorAbstraction<SecondVectorType> V2;
    typedef Common::VectorAbstraction<ThirdVectorType> V3;
    for (size_t ii = 0; ii < M::rows(QR); ++ii)
      V3::set_entry(y, ii, V2::get_entry(x, ii));
    if (false) {
      ;
#if HAVE_LAPACKE || HAVE_MKL
    } else if (has_contiguous_storage) {
      // These are the number of rows and columns of the matrix C in the documentation of dormqr.
      // As we only have a vector, i.e. C = x, the number of columns is 1.
      const size_t num_rows = x.size();
      const size_t num_cols = 1;
      auto info = Common::Lapacke::dormqr(lapacke_storage_layout(),
                                          'L',
                                          transpose == XT::Common::Transpose::yes ? 'T' : 'N',
                                          num_rows,
                                          num_cols,
                                          num_rows,
                                          M::data(QR),
                                          num_rows,
                                          V::data(tau),
                                          V3::data(y),
                                          is_row_major ? num_cols : num_rows);
      if (info)
        DUNE_THROW(Dune::MathError, "Multiplication by Q^T failed");
#endif // HAVE_LAPACKE || HAVE_MKL
    } else {
      const size_t num_rows = M::rows(QR);
      const size_t num_cols = M::cols(QR);
      VectorType w = V::create(num_rows, 0.);
      if (transpose == XT::Common::Transpose::no)
        for (int jj = num_cols - 1; jj >= 0; --jj) {
          set_w_vector(QR, jj, w);
          multiply_householder_from_left(y, tau[jj], w, jj, num_rows);
        }
      else
        for (int jj = 0; jj < int(num_cols); ++jj) {
          set_w_vector(QR, jj, w);
          multiply_householder_from_left(y, tau[jj], w, jj, num_rows);
        }
    }
  } // static void apply_q_from_qr(...)

private:
  // w is the vector [1; QR(j+1:end,j)]
  static void set_w_vector(const MatrixType& QR, const int jj, VectorType& w)
  {
    const size_t num_rows = M::rows(QR);
    V::set_entry(w, jj, 1);
    for (int ii = jj + 1; ii < int(num_rows); ++ii)
      V::set_entry(w, ii, M::get_entry(QR, ii, jj));
  }

  static typename M::template MatrixTypeTemplate<M::static_rows, M::static_rows>
  calculate_q_explicitly(const MatrixType& QR, const VectorType& tau)
  {
    const size_t num_rows = M::rows(QR);
    const size_t num_cols = M::cols(QR);
    auto ret = eye_matrix<typename M::template MatrixTypeTemplate<M::static_rows, M::static_rows>>(
        num_rows, Common::dense_pattern(num_rows, num_rows));
    VectorType w = V::create(num_rows, 1.);
    for (int jj = int(num_cols) - 1; jj >= 0; --jj) {
      set_w_vector(QR, jj, w);
      multiply_householder_from_left(ret, V::get_entry(tau, jj), w, jj, num_rows, 0, num_rows);
    }
    return ret;
  }
}; // struct QrHelper


} // namespace internal


/**
 * Performs a QR factorization with pivoting AP = QR. The matrix A is overwritten with QR, additional multipliers are
 * stored in tau and the permutations are stored in permutations.
 * \see Common::Lapacke::dgeqp3
 */
template <class MatrixType, class VectorType, class IndexVectorType>
typename std::enable_if_t<Common::is_matrix<MatrixType>::value, void>
qr(MatrixType& A, VectorType& tau, IndexVectorType& permutations)
{
  internal::QrHelper<MatrixType, VectorType, IndexVectorType>::qr(A, tau, permutations);
} // void solve_lower_triangular(...)

// calculate y =  Q * x or y = Q^T * x where Q is from the qr decomposition A = QR
template <Common::Transpose transpose,
          class MatrixType,
          class VectorType,
          class SecondVectorType,
          class ThirdVectorType>
void apply_q_from_qr(const MatrixType& QR, const VectorType& tau, const SecondVectorType& x, ThirdVectorType& y)
{
  internal::QrHelper<MatrixType, VectorType>::template apply_q_from_qr<transpose>(QR, tau, x, y);
}

// calculate y =  Q * x or y = Q^T * x where Q is from the qr decomposition A = QR
template <class MatrixType, class VectorType, class M = Common::MatrixAbstraction<MatrixType>>
typename M::template MatrixTypeTemplate<M::static_rows, M::static_rows> calculate_q_from_qr(const MatrixType& QR,
                                                                                            const VectorType& tau)
{
  return internal::QrHelper<MatrixType, VectorType>::calculate_q_from_qr(QR, tau);
}

template <class MatrixType, class IndexVectorType>
void get_permutation_matrix(const IndexVectorType& permutations, MatrixType& P)
{
  P *= 0.;
  for (size_t ii = 0; ii < permutations.size(); ++ii)
    Common::set_matrix_entry(P, permutations[ii], ii, 1.);
}

/**
 *  \brief Solves Ax = b, where AP = QR is the QR decomposition.
 *  QR has to be stored in QR and tau, P has to be stored in permutations.
 *  \see qr
 */
template <class MatrixType, class VectorType, class SecondVectorType, class RhsVectorType, class IndexVectorType>
void solve_qr_factorized(const MatrixType& QR,
                         const VectorType& tau,
                         const IndexVectorType& permutations,
                         SecondVectorType& x,
                         const RhsVectorType& b,
                         SecondVectorType* work = nullptr)
{

  typedef Common::MatrixAbstraction<MatrixType> M;
  typedef Common::VectorAbstraction<VectorType> V;
  typedef Common::VectorAbstraction<IndexVectorType> VI;

  auto num_rows = M::rows(QR);
  auto num_cols = M::cols(QR);
  if (num_rows != num_cols)
    DUNE_THROW(NotImplemented, "Not implemented for non-square matrices!");

  std::unique_ptr<SecondVectorType> work_ptr;
  if (!work) {
    work_ptr = std::make_unique<VectorType>(x);
    work = work_ptr.get();
  }
  *work = x;

  // Calculate c = Q^T b;
  apply_q_from_qr<XT::Common::Transpose::yes>(QR, tau, b, x);

  // Solve R x = c;
  if (M::storage_layout == Common::StorageLayout::csr || M::storage_layout == Common::StorageLayout::csc) {
    // Currently, the implementation of solve_upper_triangular for sparse matrices relies on a triangular pattern, so
    // copy R from QR to its own matrix first
    MatrixType R = eye_matrix<MatrixType>(
        num_rows, num_cols, Common::triangular_pattern(num_rows, num_cols, Common::MatrixPattern::upper_triangular));
    for (size_t ii = 0; ii < num_rows; ++ii)
      for (size_t jj = ii; jj < num_cols; ++jj)
        M::set_entry(R, ii, jj, M::get_entry(QR, ii, jj));
    solve_upper_triangular(R, *work, x);
  } else {
    solve_upper_triangular(QR, *work, x);
  }

  // Undo permutations
  for (int ii = 0; ii < int(M::rows(QR)); ++ii)
    V::set_entry(x, VI::get_entry(permutations, ii), V::get_entry(*work, ii));
}


} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_ALGORITHMS_QR_HH
