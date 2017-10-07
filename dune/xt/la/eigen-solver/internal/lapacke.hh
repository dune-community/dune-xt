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

#include <dune/xt/common/matrix.hh>

#include <dune/xt/la/exceptions.hh>
#include <dune/xt/la/type_traits.hh>

namespace Dune {
namespace XT {
namespace LA {
namespace internal {


#if HAVE_LAPACKE


class UnitMatrix
{
public:
  UnitMatrix(int N);
  double* get(int N);

private:
  std::unique_ptr<std::vector<double>> unit_matrix_;
}; // class UnitMatrix;

// We do not call the Lapacke functions directly to avoid including the lapacke.h header in this
// header. The lapacke header defines some macros which lead to conflicts with other includes.
struct LapackeWrapper
{
  static int dggev(char jobvl,
                   char jobvr,
                   int n,
                   double* a,
                   int lda,
                   double* b,
                   int ldb,
                   double* alphar,
                   double* alphai,
                   double* beta,
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
          bool is_matrix = XT::Common::MatrixAbstraction<MatrixType>::is_matrix>
struct lapacke_helper
{
  typedef XT::Common::MatrixAbstraction<MatrixType> MatrixAbstractionType;

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
      if (XT::Common::FloatCmp::ne(eigvals[jj].imag(), 0.)) {
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
  typedef XT::Common::MatrixAbstraction<MatrixType> MatrixAbstractionType;

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
      if (XT::Common::FloatCmp::ne(eigvals[jj].imag(), 0.)) {
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
  std::vector<double> alpha_real(N), alpha_imag(N), beta(N);

  static thread_local UnitMatrix unit_matrix(N);
  int info = LapackeWrapper::dggev('N',
                                   eigvecs ? 'V' : 'N',
                                   N,
                                   matrix,
                                   N,
                                   unit_matrix.get(N),
                                   N,
                                   alpha_real.data(),
                                   alpha_imag.data(),
                                   beta.data(),
                                   (double*)nullptr,
                                   N,
                                   eigvecs,
                                   N);

  if (info != 0)
    DUNE_THROW(Dune::MathError, "Lapack returned error " + std::to_string(info) + "!");

  for (size_t rr = 0; rr < size_t(N); ++rr) {
    assert(XT::Common::FloatCmp::ne(beta[rr], 0., 1e-6));
    eigvals[rr] = {alpha_real[rr] / beta[rr], alpha_imag[rr] / beta[rr]};
  }
} // ... compute_using_lapacke(...)

template <class Traits, class S>
std::vector<std::complex<S>> compute_all_eigenvalues_using_lapacke(const XT::LA::MatrixInterface<Traits, S>& matrix)
{
  const size_t N = matrix.rows();
  std::vector<double> tmp_matrix(N * N);
  std::vector<std::complex<S>> ret(N);
  size_t ii = 0;
  for (size_t rr = 0; rr < N; ++rr)
    for (size_t cc = 0; cc < N; ++cc)
      tmp_matrix[ii++] = matrix.get_entry(rr, cc);
  compute_using_lapacke(tmp_matrix.data(), ret, nullptr);
  return ret;
}

template <class S, int N>
std::vector<std::complex<S>> compute_all_eigenvalues_using_lapacke(const Dune::FieldMatrix<S, N, N>& matrix)
{
  std::vector<std::complex<S>> ret(N);
  std::vector<double> tmp_matrix(N * N);
  size_t ii = 0;
  for (size_t rr = 0; rr < N; ++rr)
    for (size_t cc = 0; cc < N; ++cc)
      tmp_matrix[ii++] = matrix[rr][cc];
  compute_using_lapacke(tmp_matrix.data(), ret, nullptr);
  return ret;
}

template <class Traits, class S, class ReturnType, class ReturnValueType>
void compute_all_eigenvectors_using_lapacke(const XT::LA::MatrixInterface<Traits, S>& matrix,
                                            ReturnType& ret,
                                            ReturnValueType)
{
  const size_t N = matrix.rows();
  std::vector<double> tmp_matrix(N * N), eigvecs(N * N);
  std::vector<std::complex<S>> eigvals(N);
  size_t ii = 0;
  for (size_t rr = 0; rr < N; ++rr)
    for (size_t cc = 0; cc < N; ++cc)
      tmp_matrix[ii++] = matrix.get_entry(rr, cc);
  compute_using_lapacke(tmp_matrix.data(), eigvals, eigvecs.data());
  lapacke_helper<ReturnType, S, std::is_arithmetic<ReturnValueType>::value>::set_eigvecs(ret, eigvecs.data(), eigvals);
}

template <class S, int N, class ReturnType, class ReturnValueType>
void compute_all_eigenvectors_using_lapacke(const Dune::FieldMatrix<S, N, N>& matrix, ReturnType& ret, ReturnValueType)
{
  std::vector<double> tmp_matrix(N * N), eigvecs(N * N);
  std::vector<std::complex<S>> eigvals(N);
  size_t ii = 0;
  for (size_t rr = 0; rr < N; ++rr)
    for (size_t cc = 0; cc < N; ++cc)
      tmp_matrix[ii++] = matrix[rr][cc];
  compute_using_lapacke(tmp_matrix.data(), eigvals, eigvecs.data());
  lapacke_helper<ReturnType, S, std::is_arithmetic<ReturnValueType>::value>::set_eigvecs(ret, eigvecs.data(), eigvals);
}


#else // HAVE_LAPACKE


template <class Traits, class S>
std::vector<std::complex<S>> compute_all_eigenvalues_using_lapacke(const XT::LA::MatrixInterface<Traits, S>& matrix)
{
  static_assert(AlwaysFalse<S>::value, "You are missing Lapacke!");
}

template <class S, int N>
std::vector<std::complex<S>> compute_all_eigenvalues_using_lapacke(const Dune::FieldMatrix<S, N, N>& matrix)
{
  static_assert(AlwaysFalse<S>::value, "You are missing Lapacke!");
}

template <class Traits, class S, class ReturnType, class ReturnValueType>
void compute_all_eigenvectors_using_lapacke(const XT::LA::MatrixInterface<Traits, S>& matrix,
                                            ReturnType& ret,
                                            ReturnValueType)
{
  static_assert(AlwaysFalse<S>::value, "You are missing Lapacke!");
}

template <class S, int N, class ReturnType, class ReturnValueType>
void compute_all_eigenvectors_using_lapacke(const Dune::FieldMatrix<S, N, N>& matrix, ReturnType& ret, ReturnValueType)
{
  static_assert(AlwaysFalse<S>::value, "You are missing Lapacke!");
}


#endif // HAVE_LAPACKE


} // namespace internal
} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_EIGEN_SOLVER_INTERNAL_LAPACKE_HH
