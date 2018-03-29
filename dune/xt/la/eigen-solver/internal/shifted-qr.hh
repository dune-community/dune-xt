// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017 - 2018)

#ifndef DUNE_XT_LA_EIGEN_SOLVER_INTERNAL_SHIFTEDQR_HH
#define DUNE_XT_LA_EIGEN_SOLVER_INTERNAL_SHIFTEDQR_HH

#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/common/string.hh>

#include <dune/xt/la/exceptions.hh>
#include <dune/xt/la/type_traits.hh>
#include <dune/xt/la/container/eye-matrix.hh>

namespace Dune {
namespace XT {
namespace LA {
namespace internal {


// A shifted QR eigensolver for real matrices with real eigenvalues. Complex matrices and real matrices with complex
// eigenvalues are not supported.
template <class FieldType>
struct RealQrEigenSolver
{
  typedef Dune::DynamicVector<FieldType> VectorType;
  typedef Dune::DynamicMatrix<FieldType> MatrixType;
  static constexpr size_t max_iterations = 10000;

  static std::vector<double> calculate_eigenvalues_by_shifted_qr(MatrixType& A, const std::unique_ptr<MatrixType>& Q)
  {
    const size_t num_rows = A.rows();
    const size_t num_cols = A.cols();
    const FieldType tol = 1e-15;
    auto R_k = XT::Common::make_unique<MatrixType>(num_rows, num_cols, 0.);
    auto Q_k = XT::Common::make_unique<MatrixType>(num_rows, num_cols, 0.);
    for (size_t jj = num_rows - 1; jj > 0; --jj) {
      size_t num_remaining_rows = jj + 1;
      size_t num_remaining_cols = num_remaining_rows;
      FieldType residual = std::abs(A[jj][jj - 1]);
      size_t kk = 0;
      // TODO: Choose appropiate stopping criterion
      while (XT::Common::FloatCmp::gt(residual, tol * (std::abs(A[jj][jj]) + std::abs(A[jj - 1][jj - 1])))
             && kk < max_iterations) {
        // Use Wilkinson shift, i.e. use the eigenvalue of the lower right 2x2 matrix [a b; c d] that is closer to
        // A[jj][jj].
        // If the eigenvalues are complex, we just use the eigenvalues of the symmetric matrix [a c; c d].
        // TODO: Use an appropiate shifting strategy if eigenvalues are complex
        // The eigenvalues are (a+d)/2 +- sqrt((a+d)^2/4 - (ad-bc))
        auto a = A[jj - 1][jj - 1];
        auto b = A[jj - 1][jj];
        auto c = A[jj][jj - 1];
        auto d = A[jj][jj];
        auto inside_root = (a + d) * (a + d) / 4. - (a * d - b * c);
        if (inside_root < 0.) {
          b = c;
          inside_root = (a + d) * (a + d) / 4. - (a * d - b * c);
        }
        auto eigval1 = (a + d) / 2. + std::sqrt(inside_root);
        auto eigval2 = (a + d) / 2. - std::sqrt(inside_root);
        auto shift = std::abs(eigval1 - d) < std::abs(eigval2 - d) ? eigval1 : eigval2;
        for (size_t rr = 0; rr < num_remaining_rows; ++rr)
          A[rr][rr] -= shift;

        // Calculate QR decomp. If jj < num_rows-1, A_k-shift*I has the form [A1 A2; 0 A3-shift*I],
        // so we are only calculating the QR decomposition Q1*R1 of A1.
        // Then Q_k = [Q1 0; 0 I], R_k = [R1 Q1^T*A2; 0 A3] is a QR decomp. of A_k.
        QR_decomp(A, *Q_k, *R_k, num_remaining_rows, num_remaining_cols);

        // Calculate A_{k+1} = R_k Q_k + shift I. With the QR decomposition above, this
        // is calculated as A_{k+1} = [R1*Q1+shift*I Q1^T*A2; 0 A3]
        // calculate upper left part
        for (size_t rr = 0; rr < num_remaining_rows; ++rr) {
          for (size_t cc = 0; cc < num_remaining_cols; ++cc) {
            A[rr][cc] = 0.;
            for (size_t ll = 0; ll < num_remaining_rows; ++ll)
              A[rr][cc] += (*R_k)[rr][ll] * (*Q_k)[ll][cc];
          } // cc
          A[rr][rr] += shift;
        } // rr
        // update upper right part
        // we do not need R_k anymore in this step, so use it as temporary storage
        auto& A_copy = *R_k;
        A_copy = A;
        for (size_t rr = 0; rr < num_remaining_rows; ++rr) {
          for (size_t cc = num_remaining_cols; cc < num_cols; ++cc) {
            A[rr][cc] = 0.;
            for (size_t ll = 0; ll < num_remaining_rows; ++ll)
              A[rr][cc] += (*Q_k)[ll][rr] * A_copy[ll][cc];
          } // cc
        } // rr

        if (Q) {
          // Update Q by multiplicating Q_k from the right
          auto& Q_copy = *R_k;
          Q_copy = *Q;
          for (size_t rr = 0; rr < num_rows; ++rr) {
            for (size_t cc = 0; cc < num_remaining_cols; ++cc) {
              (*Q)[rr][cc] = 0.;
              for (size_t ll = 0; ll < num_remaining_rows; ++ll)
                (*Q)[rr][cc] += Q_copy[rr][ll] * (*Q_k)[ll][cc];
            } // cc
          } // rr
        } // if (Q)

        ++kk;
        residual = std::abs(A[jj][jj - 1]);
      } // while(residual != 0)
      if (kk >= max_iterations)
        DUNE_THROW(Dune::MathError,
                   "Eigen solver did not converge (stopped after " + XT::Common::to_string(max_iterations)
                       + " iterations)");
    } // jj

    // Now eigenvalues are the diagonal elements of A
    std::vector<double> eigenvalues(num_rows);
    for (size_t rr = 0; rr < num_rows; ++rr)
      eigenvalues[rr] = A[rr][rr];
    return eigenvalues;
  } // .. calculate_eigenvalues_by_shifted_qr(...)

  // if Q is provided, Q is expected to be the unit matrix initially
  static std::vector<double> get_eigenvalues(MatrixType& A, const std::unique_ptr<MatrixType>& Q = nullptr)
  {
    hessenberg_transformation(A, Q);
    return calculate_eigenvalues_by_shifted_qr(A, Q);
  } // ... get_eigenvalues(...)

  static std::unique_ptr<MatrixType> get_eigenvectors(const MatrixType& A_triangular, const MatrixType& Q)
  {
    const MatrixType& A = A_triangular;
    auto ret = std::make_unique<MatrixType>(A);
    *ret *= 0.;
    // Calculate eigenvectors by calculating eigenvectors of triangular matrix T blockwise by backward substitution,
    // see Handbook Series Linear Algebra, Eigenvectors of Real and Complex Matrices by LR and QR triangularizations,
    // contributed by G. Peters and J. H. Wilkinsons,
    // https://link.springer.com/content/pdf/10.1007/BF02219772.pdf, equation (3)
    VectorType eigvec(A.rows());
    for (size_t kk = 0; kk < A.rows(); ++kk) {
      eigvec *= 0.; // TODO: should be unnecessary, remove
      // Compute k-th eigenvector by backsubstitution. For that purpose, set k-th component to 1 and then
      // do the usual backsubstitution. If we encounter a zero on the diagonal in the process, we just assign
      // zero to that component.
      eigvec[kk] = 1.;
      assert(kk <= std::numeric_limits<int>::max());
      for (int rr = static_cast<int>(kk) - 1; rr >= 0; --rr) {
        eigvec[rr] = 0.;
        if (XT::Common::FloatCmp::ne(A[rr][rr], A[kk][kk])) {
          for (size_t cc = rr + 1; cc <= kk; ++cc)
            eigvec[rr] -= A[rr][cc] * eigvec[cc];
          eigvec[rr] /= A[rr][rr] - A[kk][kk];
        }
      } // rr
      // apply matrix Q
      for (size_t rr = 0; rr < A.rows(); ++rr) {
        (*ret)[rr][kk] = 0.;
        for (size_t cc = 0; cc <= kk; ++cc)
          (*ret)[rr][kk] += Q[rr][cc] * eigvec[cc];
      }
    } // kk (eigenvectors)
    return ret;
  } // ... get_eigenvectors(...)

  //! \brief modified sign function returning 1 instead of 0 if the value is 0
  static FieldType xi(FieldType val)
  {
    return val < 0. ? -1. : 1.;
  }

  static FieldType get_norm_x(const MatrixType& A, const size_t col_index, size_t num_rows)
  {
    FieldType norm(0);
    for (size_t rr = col_index; rr < num_rows; ++rr)
      norm += std::pow(A[rr][col_index], 2);
    return std::sqrt(norm);
  }

  // Calculates P * A, where P = (I 0 0; 0 I-beta*u*u^T 0; 0 0 I) and u = v[first_row:past_last_row]
  static void multiply_householder_from_left(
      MatrixType& A, const FieldType& beta, const VectorType& v, const size_t first_row, const size_t past_last_row)
  {
    // calculate u^T A first
    VectorType uT_A(v.size(), 0.);
    for (size_t cc = 0; cc < A.cols(); ++cc)
      for (size_t rr = first_row; rr < past_last_row; ++rr)
        uT_A[cc] += v[rr] * A[rr][cc];
    // uT_A now contains u^T A[first_row:past_last_row,:]
    for (size_t rr = first_row; rr < past_last_row; ++rr)
      for (size_t cc = 0; cc < A.cols(); ++cc)
        A[rr][cc] -= beta * v[rr] * uT_A[cc];
  }

  // Calculates A * P.
  // \see multiply_householder_from_left
  static void multiply_householder_from_right(
      MatrixType& A, const FieldType& beta, const VectorType& v, const size_t first_col, const size_t past_last_col)
  {
    // calculate A u first
    VectorType Au(v.size(), 0.);
    for (size_t rr = 0; rr < A.rows(); ++rr)
      for (size_t cc = first_col; cc < past_last_col; ++cc)
        Au[rr] += A[rr][cc] * v[cc];
    // Au now contains A[:,first_col:past_last_col] u
    for (size_t rr = 0; rr < A.rows(); ++rr)
      for (size_t cc = first_col; cc < past_last_col; ++cc)
        A[rr][cc] -= beta * Au[rr] * v[cc];
  }

  /** \brief This is a simple QR scheme using Householder reflections.
  * The householder matrix is written as H = I - 2 v v^T, where v = u/||u|| and u = x - s ||x|| e_1, s = +-1 has the
  * opposite sign of u_1 and x is the current column of A. The matrix H is rewritten as
  * H = I - tau w w^T, where w=u/u_1 and tau = -s u_1/||x||.
  * The num_rows and num_cols parameter is used to restrict the rows and columns, the Q and R will only contain
  * the QR decomposition of A[0:num_rows, 0:num_cols].
  * \see https://en.wikipedia.org/wiki/QR_decomposition#Using_Householder_reflections.
  * \see http://www.cs.cornell.edu/~bindel/class/cs6210-f09/lec18.pdf
  * \todo The scheme does not take into account that A is in Hessenberg form, so there could be a significant
  * speedup if an according QR scheme is used (e.g. givens rotations, see https://lp.uni-goettingen.de/get/text/2138).
  */
  static void QR_decomp(const MatrixType& A, MatrixType& Q, MatrixType& R, size_t num_rows, size_t num_cols)
  {
    R = A;
    for (size_t rr = 0; rr < num_rows; ++rr) {
      Q[rr] *= 0.;
      Q[rr][rr] = 1.;
    }

    VectorType w(A.rows(), 0.);
    FieldType tau;
    for (size_t jj = 0; jj < std::min(num_rows - 1, num_cols); ++jj) {
      const auto norm_x = get_norm_x(R, jj, num_rows);
      if (XT::Common::FloatCmp::gt(norm_x, 0.)) {
        // find entry with greatest absolute value for pivoting
        size_t index = jj;
        FieldType max = std::abs(R[jj][jj]);
        for (size_t kk = jj + 1; kk < num_rows; ++kk) {
          if (XT::Common::FloatCmp::gt(std::abs(R[kk][jj]), max)) {
            max = std::abs(R[kk][jj]);
            index = kk;
          }
        }
        if (index != jj) { // swap rows
          // swapping of rows i,j can be done by Householder I - (e_i - e_j)(e_i-e_j)^T
          VectorType e_diff(A.rows(), 0.);
          e_diff[jj] = 1.;
          e_diff[index] = -1.;
          multiply_householder_from_left(R, 1., e_diff, jj, num_rows);
          multiply_householder_from_right(Q, 1., e_diff, jj, num_cols);
        }
        const auto s = -sign(R[jj][jj]);
        const FieldType u1 = R[jj][jj] - s * norm_x;
        w[jj] = 1.;
        for (size_t rr = jj + 1; rr < num_rows; ++rr)
          w[rr] = R[rr][jj] / u1;
        tau = -s * u1 / norm_x;

        // calculate R = Q_k R and Q = Q Q_k
        multiply_householder_from_left(R, tau, w, jj, num_rows);
        multiply_householder_from_right(Q, tau, w, jj, num_cols);
      } // if (norm_x != 0)
    } // jj

    // choose Q such that largest entry of each column is positive
    size_t row_index = 0;
    for (size_t cc = 0; cc < num_cols; ++cc) {
      FieldType max = std::numeric_limits<FieldType>::lowest();
      for (size_t rr = 0; rr < num_rows; ++rr) {
        if (std::abs(Q[rr][cc]) > max) {
          max = std::abs(Q[rr][cc]);
          row_index = rr;
        }
      } // rr
      if (XT::Common::FloatCmp::lt(Q[row_index][cc], 0.)) {
        // scal column of Q if largest entry is negative
        // scal row of R to ensure that still A = QR
        for (size_t rr = 0; rr < num_rows; ++rr) {
          Q[rr][cc] = -Q[rr][cc];
          R[cc][rr] = -R[cc][rr];
        } // rr
      } // if (largest entry negative)
    } // cc
  } // void QR_decomp(...)

  //! \brief Transform A to Hessenberg form by transformation P^T A P
  //! \see https://lp.uni-goettingen.de/get/text/2137
  static void hessenberg_transformation(MatrixType& A, const std::unique_ptr<MatrixType>& Q)
  {
    assert(A.rows() == A.cols() && "Hessenberg transformation needs a square matrix!");
    VectorType u(A.rows(), 0.);
    for (size_t jj = 0; jj < A.rows() - 2; ++jj) {
      FieldType gamma = 0;
      for (size_t rr = jj + 1; rr < A.rows(); ++rr)
        gamma += std::pow(A[rr][jj], 2);
      gamma = std::sqrt(gamma);
      FieldType beta = gamma * (gamma + std::abs(A[jj + 1][jj]));
      if (XT::Common::FloatCmp::ne(gamma, 0.) && XT::Common::FloatCmp::ne(beta, 0.)) {
        beta = 1. / beta;
        for (size_t rr = jj + 1; rr < A.rows(); ++rr)
          u[rr] = A[rr][jj];
        u[jj + 1] += xi(A[jj + 1][jj]) * gamma;
        // calculate P A P with P = diag(I_j, (I_{n-j} - beta u u*))
        // calculate P A = A - (beta u) (u^T A) first
        multiply_householder_from_left(A, beta, u, jj + 1, A.rows());
        // now calculate (PA) P  = PA - (PA u) (beta u^T)
        multiply_householder_from_right(A, beta, u, jj + 1, A.cols());
        if (Q)
          multiply_householder_from_right(*Q, beta, u, jj + 1, A.cols());
      } // if (gamma != 0)
    } // jj
  } // void hessenberg_transformation(...)
}; // class RealQrEigenSolver<...>

template <class FieldType>
constexpr size_t RealQrEigenSolver<FieldType>::max_iterations;

template <class MatrixType>
Dune::DynamicMatrix<typename Common::MatrixAbstraction<MatrixType>::RealType>
copy_to_dynamic_matrix(const MatrixType& matrix)
{
  typedef XT::Common::MatrixAbstraction<MatrixType> MatAbstrType;
  // copy matrix to DynamicMatrix
  size_t num_rows = MatAbstrType::rows(matrix);
  size_t num_cols = MatAbstrType::cols(matrix);
  Dune::DynamicMatrix<typename MatAbstrType::RealType> ret(num_rows, num_cols, 0.);
  for (size_t ii = 0; ii < num_rows; ++ii)
    for (size_t jj = 0; jj < num_cols; ++jj)
      ret[ii][jj] = MatAbstrType::get_entry(matrix, ii, jj);
  return ret;
}

template <class MatrixType>
void copy_from_dynamic_matrix(
    const Dune::DynamicMatrix<typename Common::MatrixAbstraction<MatrixType>::RealType>& dyn_matrix, MatrixType& matrix)
{
  typedef XT::Common::MatrixAbstraction<MatrixType> MatAbstrType;
  size_t num_rows = MatAbstrType::rows(matrix);
  size_t num_cols = MatAbstrType::cols(matrix);
  assert(dyn_matrix.rows() == num_rows && dyn_matrix.cols() == num_cols);
  for (size_t ii = 0; ii < num_rows; ++ii)
    for (size_t jj = 0; jj < num_cols; ++jj)
      MatAbstrType::set_entry(matrix, ii, jj, dyn_matrix[ii][jj]);
}

template <class MatrixType>
typename std::enable_if<Common::is_matrix<MatrixType>::value, std::vector<double>>::type
compute_eigenvalues_using_qr(const MatrixType& matrix)
{
  auto tmp_matrix = copy_to_dynamic_matrix(matrix);
  return RealQrEigenSolver<typename XT::Common::MatrixAbstraction<MatrixType>::RealType>::get_eigenvalues(tmp_matrix);
}

template <class MatrixType, class EigenVectorType>
typename std::enable_if<Common::is_matrix<MatrixType>::value && Common::is_matrix<EigenVectorType>::value, void>::type
compute_real_eigenvalues_and_real_right_eigenvectors_using_qr(const MatrixType& matrix,
                                                              std::vector<double>& eigenvalues,
                                                              EigenVectorType& right_eigenvectors)
{
  auto tmp_matrix = copy_to_dynamic_matrix(matrix);
  typedef typename Dune::DynamicMatrix<typename Common::MatrixAbstraction<MatrixType>::RealType> DynamicMatrixType;
  auto Q =
      std::make_unique<DynamicMatrixType>(XT::LA::eye_matrix<DynamicMatrixType>(tmp_matrix.rows(), tmp_matrix.cols()));
  eigenvalues =
      RealQrEigenSolver<typename XT::Common::MatrixAbstraction<MatrixType>::RealType>::get_eigenvalues(tmp_matrix, Q);
  auto tmp_eigenvectors =
      RealQrEigenSolver<typename XT::Common::MatrixAbstraction<MatrixType>::RealType>::get_eigenvectors(tmp_matrix, *Q);
  copy_from_dynamic_matrix(*tmp_eigenvectors, right_eigenvectors);
}


} // namespace internal
} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_EIGEN_SOLVER_INTERNAL_SHIFTEDQR_HH
