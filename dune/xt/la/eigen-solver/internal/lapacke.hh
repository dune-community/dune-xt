// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2017 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler  (2017)

#ifndef DUNE_XT_LA_EIGEN_SOLVER_INTERNAL_LAPACKE_HH
#define DUNE_XT_LA_EIGEN_SOLVER_INTERNAL_LAPACKE_HH

#include <string>

#include <dune/common/typetraits.hh>

#include <dune/xt/common/matrix.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/la/exceptions.hh>
#include <dune/xt/la/container/matrix-interface.hh>
#include <dune/xt/la/type_traits.hh>

namespace Dune {
namespace XT {
namespace LA {
namespace internal {

#if HAVE_LAPACKE


// We do not call the Lapacke functions directly to avoid including the lapacke.h header in this
// header. The lapacke header defines some macros which lead to conflicts with other includes.
struct LapackeWrapper
{
  static int dgeev(char jobvl,
                   char jobvr,
                   int n,
                   double* a,
                   int lda,
                   double* wr,
                   double* wi,
                   double* vl,
                   int ldvl,
                   double* vr,
                   int ldvr);
};


// if MatrixType is actually not a matrix but a vector of column vectors, the indices are the other way around (columns
// first, then rows)
template <class MatrixType,
          class S,
          bool real_eigvecs_requested,
          bool is_matrix = Common::MatrixAbstraction<MatrixType>::is_matrix>
struct lapacke_helper
{
  typedef Common::MatrixAbstraction<MatrixType> MatrixAbstractionType;

  static void set_eigvecs(MatrixType& eigvecs, double* eigvecs_double, std::vector<std::complex<S>>& eigvals)
  {
    // from lapacke documentation:
    // the right eigenvectors v(j) are stored one after another in the columns of VR,
    // in the same order as their eigenvalues. If the j-th eigenvalue is real, then
    // v(j) = VR(:,j), the j-th column of VR. If the j-th and j+1)-th eigenvalues form
    // a complex conjugate pair, then v(j) = VR(:,j)+i*VR(:,j+1) and
    // v(j+1) = VR(:,j)-i*VR(:,j+1).
    const auto N = eigvals.size();
    for (size_t jj = 0; jj < N; ++jj)
      if (Common::FloatCmp::ne(eigvals[jj].imag(), 0.)) {
        for (size_t ii = 0; ii < N; ++ii) {
          MatrixAbstractionType::set_entry(
              eigvecs, ii, jj, {*(eigvecs_double + (N * ii + jj)), *(eigvecs_double + (N * ii + jj + 1))});
          MatrixAbstractionType::set_entry(
              eigvecs, ii, jj + 1, {*(eigvecs_double + (N * ii + jj)), -*(eigvecs_double + (N * ii + jj + 1))});
        }
        ++jj;
      } else {
        for (size_t ii = 0; ii < N; ++ii)
          MatrixAbstractionType::set_entry(eigvecs, ii, jj, {*(eigvecs_double + (N * ii + jj)), 0.});
      }
  }
};

template <class MatrixType, class S>
struct lapacke_helper<MatrixType, S, true, true>
{
  typedef Common::MatrixAbstraction<MatrixType> MatrixAbstractionType;

  static void set_eigvecs(MatrixType& eigvecs, double* eigvecs_double, std::vector<std::complex<S>>& eigvals)
  {
    const auto N = eigvals.size();
    for (size_t ii = 0; ii < N; ++ii)
      for (size_t jj = 0; jj < N; ++jj)
        MatrixAbstractionType::set_entry(eigvecs, ii, jj, *(eigvecs_double + (N * ii + jj)));
  }
};

template <class MatrixType, class S>
struct lapacke_helper<MatrixType, S, false, false>
{
  static void set_eigvecs(MatrixType& eigvecs, double* eigvecs_double, std::vector<std::complex<S>>& eigvals)
  {
    // from lapacke documentation:
    // the right eigenvectors v(j) are stored one after another in the columns of VR,
    // in the same order as their eigenvalues. If the j-th eigenvalue is real, then
    // v(j) = VR(:,j), the j-th column of VR. If the j-th and j+1)-th eigenvalues form
    // a complex conjugate pair, then v(j) = VR(:,j)+i*VR(:,j+1) and
    // v(j+1) = VR(:,j)-i*VR(:,j+1).
    const auto N = eigvals.size();
    for (size_t jj = 0; jj < N; ++jj)
      if (Common::FloatCmp::ne(eigvals[jj].imag(), 0.)) {
        for (size_t ii = 0; ii < N; ++ii) {
          eigvecs[jj][ii] = {*(eigvecs_double + (N * ii + jj)), *(eigvecs_double + (N * ii + jj + 1))};
          eigvecs[jj + 1][ii] = {*(eigvecs_double + (N * ii + jj)), -*(eigvecs_double + (N * ii + jj + 1))};
        }
        ++jj;
      } else {
        for (size_t ii = 0; ii < N; ++ii)
          eigvecs[jj][ii] = {*(eigvecs_double + (N * ii + jj)), 0.};
      }
  }
};

template <class MatrixType, class S>
struct lapacke_helper<MatrixType, S, true, false>
{
  static void set_eigvecs(MatrixType& eigvecs, double* eigvecs_double, std::vector<std::complex<S>>& eigvals)
  {
    const auto N = eigvals.size();
    for (size_t ii = 0; ii < N; ++ii)
      for (size_t jj = 0; jj < N; ++jj)
        eigvecs[jj][ii] = *(eigvecs_double + (N * ii + jj));
  }
};


template <class S>
void compute_using_lapacke(double* matrix, std::vector<std::complex<S>>& eigvals, double* eigvecs)
{
  int N = (int)eigvals.size();
  std::vector<double> wr(N), wi(N);

  int info = LapackeWrapper::dgeev(
      'N', eigvecs ? 'V' : 'N', N, matrix, N, wr.data(), wi.data(), (double*)nullptr, N, eigvecs, N);

  if (info != 0)
    DUNE_THROW(Exceptions::eigen_solver_failed, "Lapack returned error " + std::to_string(info) + "!");

  for (size_t rr = 0; rr < size_t(N); ++rr)
    eigvals[rr] = {wr[rr], wi[rr]};
} // ... compute_using_lapacke(...)


template <class Traits, class S>
typename std::enable_if<Common::is_arithmetic<S>::value, std::vector<Common::complex_t<S>>>::type
compute_all_eigenvalues_using_lapacke(const XT::LA::MatrixInterface<Traits, S>& matrix)
{
  const size_t N = matrix.rows();
  if (matrix.cols() != N)
    DUNE_THROW(Exceptions::eigen_solver_failed_bc_data_did_not_fulfill_requirements,
               "Matrix has to be square, is " << N << "x" << matrix.cols() << "!");
  std::vector<Common::complex_t<S>> ret(N);
  compute_using_lapacke(Common::serialize_rowwise<double>(matrix.as_imp()).get(), ret, nullptr);
  return ret;
} // ... compute_all_eigenvalues_using_lapacke(...)

template <class M>
typename std::enable_if<Common::is_matrix<M>::value && !is_matrix<M>::value,
                        std::vector<Common::complex_t<typename Common::MatrixAbstraction<M>::RealType>>>::type
compute_all_eigenvalues_using_lapacke(const M& matrix)
{
  const size_t N = Common::get_matrix_rows(matrix);
  if (Common::get_matrix_cols(matrix) != N)
    DUNE_THROW(Exceptions::eigen_solver_failed_bc_data_did_not_fulfill_requirements,
               "Matrix has to be square, is " << N << "x" << Common::get_matrix_cols(matrix) << "!");
  std::vector<Common::complex_t<typename Common::MatrixAbstraction<M>::RealType>> ret(N);
  compute_using_lapacke(Common::serialize_rowwise<double>(matrix).get(), ret, nullptr);
  return ret;
} // ... compute_all_eigenvalues_using_lapacke(...)


template <class Traits, class S, class EV_Traits>
void compute_all_eigenvectors_using_lapacke(const XT::LA::MatrixInterface<Traits, S>& matrix,
                                            XT::LA::MatrixInterface<EV_Traits, Common::complex_t<S>>& eigenvectors)
{
  const size_t N = matrix.rows();
  if (matrix.cols() != N)
    DUNE_THROW(Exceptions::eigen_solver_failed_bc_data_did_not_fulfill_requirements,
               "Matrix has to be square, is " << N << "x" << matrix.cols() << "!");
  std::vector<Common::complex_t<S>> tmp_eigenvalues(N);
  std::vector<double> tmp_eigenvectors(N * N);
  compute_using_lapacke(
      Common::serialize_rowwise<double>(matrix.as_imp()).get(), tmp_eigenvalues, tmp_eigenvectors.data());
  lapacke_helper<typename XT::LA::MatrixInterface<EV_Traits, Common::complex_t<S>>::derived_type, S, false>::
      set_eigvecs(eigenvectors.as_imp(), tmp_eigenvectors.data(), tmp_eigenvalues);
}

template <class MatrixType, class EigenValueMatrixType>
typename std::enable_if<Common::is_matrix<MatrixType>::value && !is_matrix<MatrixType>::value
                            && (Common::is_matrix<EigenValueMatrixType>::value
                                && !is_matrix<EigenValueMatrixType>::value),
                        void>::type
compute_all_eigenvectors_using_lapacke(const MatrixType& matrix, EigenValueMatrixType& eigenvectors)
{
  const size_t N = Common::get_matrix_rows(matrix);
  if (Common::get_matrix_cols(matrix) != N)
    DUNE_THROW(Exceptions::eigen_solver_failed_bc_data_did_not_fulfill_requirements,
               "Matrix has to be square, is " << N << "x" << Common::get_matrix_cols(matrix) << "!");
  if (Common::get_matrix_rows(eigenvectors) != N || Common::get_matrix_cols(eigenvectors) != N)
    DUNE_THROW(Exceptions::eigen_solver_failed_bc_data_did_not_fulfill_requirements,
               "Eigenvalue matrix should be " << N << "x" << Common::get_matrix_cols(eigenvectors) << ", is "
                                              << Common::get_matrix_rows(eigenvectors)
                                              << "x"
                                              << Common::get_matrix_cols(eigenvectors)
                                              << "!");
  using S = typename Common::MatrixAbstraction<EigenValueMatrixType>::S;
  std::vector<std::complex<Common::real_t<S>>> tmp_eigenvalues(N);
  std::vector<double> tmp_eigenvectors(N * N);
  compute_using_lapacke(Common::serialize_rowwise<double>(matrix).get(), tmp_eigenvalues, tmp_eigenvectors.data());
  lapacke_helper<EigenValueMatrixType, Common::real_t<S>, Common::is_arithmetic<S>::value>::set_eigvecs(
      eigenvectors, tmp_eigenvectors.data(), tmp_eigenvalues);
} // ... compute_all_eigenvalues_using_lapacke(...)


// template <class S, int N, class ReturnType, class ReturnValueType>
// void compute_all_eigenvectors_using_lapacke(const Dune::FieldMatrix<S, N, N>& matrix, ReturnType& ret,
// ReturnValueType)
//{
//  std::vector<double> tmp_matrix(N * N), eigvecs(N * N);
//  std::vector<std::complex<S>> eigvals(N);
//  size_t ii = 0;
//  for (size_t rr = 0; rr < N; ++rr)
//    for (size_t cc = 0; cc < N; ++cc)
//      tmp_matrix[ii++] = matrix[rr][cc];
//  compute_using_lapacke(tmp_matrix.data(), eigvals, eigvecs.data());
//  lapacke_helper<ReturnType, S, std::is_arithmetic<ReturnValueType>::value>::set_eigvecs(ret, eigvecs.data(),
//  eigvals);
//}


#else // HAVE_LAPACKE


template <class Traits, class S>
typename std::enable_if<Common::is_arithmetic<S>::value, std::vector<Common::complex_t<S>>>::type
compute_all_eigenvalues_using_lapacke(const XT::LA::MatrixInterface<Traits, S>& /*matrix*/)
{
  static_assert(AlwaysFalse<Traits>::value, "You are missing lapacke!");
}


template <class M>
typename std::enable_if<Common::is_matrix<M>::value && !is_matrix<M>::value,
                        std::vector<Common::complex_t<typename Common::MatrixAbstraction<M>::RealType>>>::type
compute_all_eigenvalues_using_lapacke(const M& /*matrix*/)
{
  static_assert(AlwaysFalse<M>::value, "You are missing lapacke!");
}

template <class Traits, class S, class EV_Traits>
void compute_all_eigenvectors_using_lapacke(const XT::LA::MatrixInterface<Traits, S>& /*matrix*/,
                                            XT::LA::MatrixInterface<EV_Traits, Common::complex_t<S>>& /*eigenvectors*/)
{
  static_assert(AlwaysFalse<Traits>::value, "You are missing lapacke!");
}

template <class MatrixType, class EigenValueMatrixType>
typename std::enable_if<Common::is_matrix<MatrixType>::value && !is_matrix<MatrixType>::value
                            && (Common::is_matrix<EigenValueMatrixType>::value
                                && !is_matrix<EigenValueMatrixType>::value),
                        void>::type
compute_all_eigenvectors_using_lapacke(const MatrixType& /*matrix*/, EigenValueMatrixType& /*eigenvectors*/)
{
  static_assert(AlwaysFalse<MatrixType>::value, "You are missing lapacke!");
}


#endif // HAVE_LAPACKE


} // namespace internal
} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_EIGEN_SOLVER_INTERNAL_LAPACKE_HH
