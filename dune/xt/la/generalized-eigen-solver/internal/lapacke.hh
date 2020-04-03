// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)
//   René Fritze     (2019)
//   Tobias Leibner  (2020)

#ifndef DUNE_XT_LA_GENERALIZED_EIGEN_SOLVER_INTERNAL_LAPACKE_HH
#define DUNE_XT_LA_GENERALIZED_EIGEN_SOLVER_INTERNAL_LAPACKE_HH

#include <complex>
#include <vector>
#include <string>
#include <numeric>
#include <type_traits>

#include <dune/common/typetraits.hh>

#include <dune/xt/common/matrix.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/common/lapacke.hh>

#include <dune/xt/la/container/common/matrix/dense.hh>
#include <dune/xt/la/exceptions.hh>
#include <dune/xt/la/type_traits.hh>

#include <dune/xt/la/eigen-solver/internal/lapacke.hh>

namespace Dune {
namespace XT {
namespace LA {
namespace internal {


/**
 * \sa https://software.intel.com/en-us/mkl-developer-reference-c-sygv
 * \note Most likely, you do not want to use this function directly, but compute_generalized_eigenvalues_using_lapack.
 */
template <class RealMatrixType>
typename std::enable_if<Common::is_matrix<std::decay_t<RealMatrixType>>::value, std::vector<std::complex<double>>>::type
compute_generalized_eigenvalues_of_real_matrices_using_lapack(RealMatrixType&& lhs_matrix, RealMatrixType&& rhs_matrix)
{
  MatrixDataProvider<std::decay_t<RealMatrixType>, is_contiguous_and_mutable<RealMatrixType>::value>
      lhs_matrix_data_provider(lhs_matrix);
  MatrixDataProvider<std::decay_t<RealMatrixType>, is_contiguous_and_mutable<RealMatrixType>::value>
      rhs_matrix_data_provider(rhs_matrix);
  return compute_generalized_eigenvalues_of_real_matrices_using_lapack_impl(
      lhs_matrix, lhs_matrix_data_provider, rhs_matrix, rhs_matrix_data_provider);
}


template <class RealMatrixType, bool contiguous_and_mutable>
typename std::enable_if<Common::is_matrix<RealMatrixType>::value, std::vector<std::complex<double>>>::type
compute_generalized_eigenvalues_of_real_matrices_using_lapack_impl(
    const RealMatrixType& lhs_matrix,
    MatrixDataProvider<RealMatrixType, contiguous_and_mutable>& lhs_matrix_data_provider,
    const RealMatrixType& rhs_matrix,
    MatrixDataProvider<RealMatrixType, contiguous_and_mutable>& rhs_matrix_data_provider)
{
  if (!Common::Lapacke::available())
    DUNE_THROW(Exceptions::generalized_eigen_solver_failed_bc_it_was_not_set_up_correctly,
               "Do not call any lapack related method if Common::Lapacke::available() is false!");
  using real_type = typename Dune::XT::Common::MatrixAbstraction<RealMatrixType>::S;
  static_assert(Dune::XT::Common::is_arithmetic<real_type>::value && !Dune::XT::Common::is_complex<real_type>::value,
                "Not implemented for complex matrices (yet)!");
  const size_t size = Dune::XT::Common::get_matrix_rows(lhs_matrix);
  const auto sz = XT::Common::numeric_cast<int>(size);
#ifdef DUNE_XT_LA_DISABLE_ALL_CHECKS
  assert(Dune::XT::Common::get_matrix_cols(lhs_matrix) == size);
  assert(Dune::XT::Common::get_matrix_rows(rhs_matrix) == size);
  assert(Dune::XT::Common::get_matrix_cols(rhs_matrix) == size);
#else
  if (Dune::XT::Common::get_matrix_cols(lhs_matrix) != size)
    DUNE_THROW(Exceptions::generalized_eigen_solver_failed_bc_data_did_not_fulfill_requirements,
               "Given LHS matrix has to be square, is " << size << "x" << Dune::XT::Common::get_matrix_cols(lhs_matrix)
                                                        << "!");
  if (Dune::XT::Common::get_matrix_rows(rhs_matrix) != size)
    DUNE_THROW(Exceptions::generalized_eigen_solver_failed_bc_data_did_not_fulfill_requirements,
               "Given matrices have to be of same size, are " << size << "x"
                                                              << Dune::XT::Common::get_matrix_cols(lhs_matrix) << "and "
                                                              << Dune::XT::Common::get_matrix_rows(rhs_matrix) << "x"
                                                              << Dune::XT::Common::get_matrix_cols(rhs_matrix) << "!");
  if (Dune::XT::Common::get_matrix_cols(rhs_matrix) != size)
    DUNE_THROW(Exceptions::generalized_eigen_solver_failed_bc_data_did_not_fulfill_requirements,
               "Given RHS matrix has to be square, is " << Dune::XT::Common::get_matrix_rows(rhs_matrix) << "x"
                                                        << Dune::XT::Common::get_matrix_cols(rhs_matrix) << "!");
#endif // DUNE_XT_LA_DISABLE_ALL_CHECKS
  thread_local std::vector<double> real_part_of_eigenvalues(size, 0.);
  int storage_layout = MatrixDataProvider<RealMatrixType, contiguous_and_mutable>::storage_layout
                               == Common::StorageLayout::dense_row_major
                           ? Common::Lapacke::row_major()
                           : Common::Lapacke::col_major();
  const int info = XT::Common::Lapacke::dsygv(storage_layout,
                                              /*problem type=*/1,
                                              /*only eigenvalues*/ 'N',
                                              /*upper triangles*/ 'U',
                                              sz,
                                              lhs_matrix_data_provider.data(),
                                              sz,
                                              rhs_matrix_data_provider.data(),
                                              sz,
                                              real_part_of_eigenvalues.data());
  if (info != 0)
    DUNE_THROW(Dune::XT::LA::Exceptions::generalized_eigen_solver_failed,
               "The lapack backend reported '"
                   << info << "', see https://software.intel.com/en-us/mkl-developer-reference-c-sygv!");
  std::vector<std::complex<double>> eigenvalues(size);
  for (size_t ii = 0; ii < size; ++ii)
    eigenvalues[ii] = {real_part_of_eigenvalues[ii], 0.};
  return eigenvalues;
} // ... compute_generalized_eigenvalues_of_real_matrices_using_lapack_impl(...)


/**
 * \note Most likely, you do not want to use this function directly, but compute_generalized_eigenvalues_using_lapack.
 */
template <class MatrixType>
struct generalized_eigenvalues_lapack_helper
{
  static_assert(Common::is_matrix<MatrixType>::value, "");

  template <bool is_complex = Common::is_complex<typename Common::MatrixAbstraction<MatrixType>::S>::value,
            bool anything = true>
  struct dtype_switch;

  template <bool anything>
  struct dtype_switch<true, anything>
  {
    template <class MatrixImp>
    static inline std::vector<std::complex<double>> eigenvalues(MatrixImp&& /*lhs_matrix*/, MatrixImp&& /*rhs_matrix*/)
    {
      static_assert(AlwaysFalse<MatrixImp>::value,
                    "Not yet implemented for complex matrices, take a look at "
                    "https://software.intel.com/en-us/mkl-developer-reference-c-sygv "
                    "and add a corresponding free function like "
                    "compute_generalized_eigenvalues_of_real_matrices_using_lapack(...)!");
      return std::vector<std::complex<double>>();
    }
  };

  template <bool anything>
  struct dtype_switch<false, anything>
  {
    template <class MatrixImp>
    static inline std::vector<std::complex<double>> eigenvalues(MatrixImp&& lhs_matrix, MatrixImp&& rhs_matrix)
    {
      return compute_generalized_eigenvalues_of_real_matrices_using_lapack(std::forward<MatrixImp>(lhs_matrix),
                                                                           std::forward<MatrixImp>(rhs_matrix));
    }
  };
}; // class generalized_eigenvalues_lapack_helper


template <class MatrixType>
typename std::enable_if<Common::is_matrix<std::decay_t<MatrixType>>::value, std::vector<std::complex<double>>>::type
compute_generalized_eigenvalues_using_lapack(MatrixType&& lhs_matrix, MatrixType&& rhs_matrix)
{
  return generalized_eigenvalues_lapack_helper<std::decay_t<MatrixType>>::template dtype_switch<>::eigenvalues(
      std::forward<MatrixType>(lhs_matrix), std::forward<MatrixType>(rhs_matrix));
}


} // namespace internal
} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_GENERALIZED_EIGEN_SOLVER_INTERNAL_LAPACKE_HH
