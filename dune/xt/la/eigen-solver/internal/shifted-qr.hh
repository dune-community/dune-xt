// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2017 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner  (2017-2018)

#ifndef DUNE_XT_LA_EIGEN_SOLVER_INTERNAL_SHIFTEDQR_HH
#define DUNE_XT_LA_EIGEN_SOLVER_INTERNAL_SHIFTEDQR_HH

#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/common/string.hh>

#include <dune/xt/la/exceptions.hh>
#include <dune/xt/la/type_traits.hh>

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

  static std::vector<double> calculate_eigenvalues_by_shifted_qr(MatrixType& A)
  {
    const size_t num_rows = A.rows();
    const size_t num_cols = A.cols();
    static constexpr size_t max_iterations = 10000;
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

        // calculate QR decomp. If jj < num_rows-1, A has the form [A1 A2; 0 A3],
        // so we are only calculating the QR decomposition Q1*R1 of A1.
        // Then Q = [Q1 0; 0 I], R = [R1 Q1^T*A2; 0 A3] is a QR decomp. of A.
        QR_decomp(A, *Q_k, *R_k, num_remaining_rows, num_remaining_cols);

        // calculate A_{k+1} = R_k Q_k + shift I. We are only interested in the diagonal
        // elements of A, and we do not reuse the other parts of A, so we are only
        // updating the upper left part.
        for (size_t rr = 0; rr < num_remaining_rows; ++rr) {
          for (size_t cc = 0; cc < num_remaining_cols; ++cc) {
            A[rr][cc] = 0.;
            for (size_t ll = 0; ll < num_remaining_rows; ++ll)
              A[rr][cc] += (*R_k)[rr][ll] * (*Q_k)[ll][cc];
          } // cc
        } // rr

        // we do not need R_k anymore in this step, so use it as temporary storage
        auto& A_copy = *R_k;
        A_copy = A;
        // update upper right part
        for (size_t rr = 0; rr < num_remaining_rows; ++rr) {
          for (size_t cc = num_remaining_cols; cc < num_cols; ++cc) {
            A[rr][cc] = 0.;
            for (size_t ll = 0; ll < num_remaining_rows; ++ll)
              A[rr][cc] += (*Q_k)[ll][rr] * A_copy[ll][cc];
          } // cc
        } // rr

        for (size_t rr = 0; rr < num_remaining_rows; ++rr)
          A[rr][rr] += shift;
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

  static void reduced_row_echelon_form(MatrixType& A, std::vector<size_t>& pivot_variables)
  {
    pivot_variables.clear();
    size_t col = 0;
    for (size_t row = 0; row < A.rows(); ++row) {
      while (col < A.cols()) {
        // find pivot
        size_t pivot = row;
        for (size_t rr = row + 1; rr < A.rows(); ++rr)
          if (std::abs(A[rr][col]) > std::abs(A[pivot][col]))
            pivot = rr;
        // if all entries in column are zero, continue with next column
        if (XT::Common::FloatCmp::eq(A[pivot][col], 0.)) {
          ++col;
          continue;
        } else {
          if (pivot != row) // swap rows
            for (size_t cc = col; cc < A.cols(); ++cc)
              std::swap(A[row][cc], A[pivot][cc]);
          // divide row by pivot element to make pivot element 1
          for (size_t cc = col + 1; cc < A.cols(); ++cc)
            A[row][cc] /= A[row][col];
          A[row][col] = 1.;
          // add row to all other rows to zero out other entries in column
          for (size_t rr = 0; rr < A.rows(); ++rr) {
            auto factor = A[rr][col];
            if (rr != row && XT::Common::FloatCmp::ne(factor, 0.)) {
              for (size_t cc = col + 1; cc < A.cols(); ++cc)
                A[rr][cc] -= A[row][cc] * factor;
              A[rr][col] = 0.;
            }
          } // rr
          // store column as pivot variable
          pivot_variables.push_back(col);
          // continue with next row and col
          ++col;
          break;
        }
      } // col
    } // row
  } // void reduced_row_echelon_form(...)

  static std::vector<double> get_eigenvalues(const MatrixType& A_in)
  {
    auto A_ptr = std::make_unique<MatrixType>(A_in);
    auto& A = *A_ptr;
    hessenberg_transformation(A);
    auto eigenvalues = calculate_eigenvalues_by_shifted_qr(A);
    std::sort(eigenvalues.begin(), eigenvalues.end());
    return eigenvalues;
  } // ... get_eigenvalues(...)

  static std::unique_ptr<MatrixType> get_eigenvectors(const MatrixType& A_in, const std::vector<double>& eigenvalues)
  {
    auto A_ptr = std::make_unique<MatrixType>(A_in);
    auto& A = *A_ptr;
    auto ret = std::make_unique<MatrixType>(A_in);
    // Calculate eigenvectors by solving (A - \lambda I) x = 0.
    for (size_t ii = 0; ii < eigenvalues.size(); ++ii) {
      // get eigenvalue, calculate multiplicity
      double lambda = eigenvalues[ii];
      size_t multiplicity = 1;
      while (XT::Common::FloatCmp::eq(eigenvalues[ii], eigenvalues[ii + 1])) {
        ++multiplicity;
        ++ii;
      }
      // get matrix B = A - \lambda I
      A = A_in;
      for (size_t rr = 0; rr < A.rows(); ++rr)
        A[rr][rr] -= lambda;

      // transform B to reduced row echelon form
      // in the process, keep track of pivot variables
      std::vector<size_t> pivot_variables;
      reduced_row_echelon_form(A, pivot_variables);

      // Now get the eigenvectors from the reduced row echelon form. To get an eigenvector, set one of the free
      // variables to 1 and the other free variables to 0 and calculate the pivot variables.
      std::vector<size_t> free_variables;
      for (size_t jj = 0; jj < A.cols(); ++jj)
        if (std::find(pivot_variables.begin(), pivot_variables.end(), jj) == pivot_variables.end())
          free_variables.push_back(jj);
      if (free_variables.size() != multiplicity)
        DUNE_THROW(Dune::MathError, "Failed to compute eigenvectors!");
      size_t eigenvector_index = ii - multiplicity + 1;
      for (const auto& cc : free_variables) {
        for (size_t rr = 0; rr < A.rows(); ++rr)
          (*ret)[rr][eigenvector_index] = -A[rr][cc];
        (*ret)[cc][eigenvector_index] = 1.;
        ++eigenvector_index;
      } // cc
    } // ii
    // normalize eigenvectors
    for (size_t col = 0; col < ret->cols(); ++col) {
      double norm = 0.;
      for (size_t row = 0; row < ret->rows(); ++row)
        norm += std::pow((*ret)[row][col], 2);
      auto inv_norm = 1. / std::sqrt(norm);
      for (size_t row = 0; row < ret->rows(); ++row)
        (*ret)[row][col] *= inv_norm;
    }
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
  static void hessenberg_transformation(MatrixType& A)
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
      } // if (gamma != 0)
    } // jj
  } // void hessenberg_transformation(...)
}; // class RealQrEigenSolver<...>

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
  return RealQrEigenSolver<typename MatrixType::FieldType>::get_eigenvalues(tmp_matrix);
}

template <class MatrixType, class EigenVectorType>
typename std::enable_if<Common::is_matrix<MatrixType>::value && Common::is_matrix<EigenVectorType>::value, void>::type
compute_real_eigenvalues_and_real_right_eigenvectors_using_qr(const MatrixType& matrix,
                                                              std::vector<double>& eigenvalues,
                                                              EigenVectorType& right_eigenvectors)
{
  auto tmp_matrix = copy_to_dynamic_matrix(matrix);
  eigenvalues =
      RealQrEigenSolver<typename XT::Common::MatrixAbstraction<MatrixType>::RealType>::get_eigenvalues(tmp_matrix);
  auto tmp_eigenvectors =
      RealQrEigenSolver<typename XT::Common::MatrixAbstraction<MatrixType>::RealType>::get_eigenvectors(tmp_matrix,
                                                                                                        eigenvalues);
  copy_from_dynamic_matrix(*tmp_eigenvectors, right_eigenvectors);
}


} // namespace internal
} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_EIGEN_SOLVER_INTERNAL_SHIFTEDQR_HH
