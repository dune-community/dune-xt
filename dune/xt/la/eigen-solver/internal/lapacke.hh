// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)
//   Tobias Leibner  (2017 - 2018)

#ifndef DUNE_XT_LA_EIGEN_SOLVER_INTERNAL_LAPACKE_HH
#define DUNE_XT_LA_EIGEN_SOLVER_INTERNAL_LAPACKE_HH

#include <complex>
#include <vector>
#include <string>
#include <numeric>

#include <dune/common/typetraits.hh>

#include <dune/xt/common/matrix.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/common/lapacke.hh>
#include <dune/xt/la/exceptions.hh>
#include <dune/xt/la/container/matrix-interface.hh>
#include <dune/xt/la/type_traits.hh>

namespace Dune {
namespace XT {
namespace LA {
namespace internal {


/**
 * \sa https://software.intel.com/en-us/mkl-developer-reference-c-geev
 * \note Most likely, you do not want to use this function directly, but compute_eigenvalues_using_lapack.
 */
template <class SerializableRealMatrixType>
typename std::enable_if<Common::is_matrix<SerializableRealMatrixType>::value, std::vector<std::complex<double>>>::type
compute_eigenvalues_of_a_real_matrix_using_lapack(const SerializableRealMatrixType& serializable_matrix)
{
  if (!Common::Lapacke::available())
    DUNE_THROW(Exceptions::eigen_solver_failed_bc_it_was_not_set_up_correctly,
               "Do not call any lapack related method if Common::Lapacke::available() is false!");
  using real_type = typename Dune::XT::Common::MatrixAbstraction<SerializableRealMatrixType>::S;
  static_assert(Dune::XT::Common::is_arithmetic<real_type>::value && !Dune::XT::Common::is_complex<real_type>::value,
                "Not implemented for complex matrices (yet)!");
  const size_t size = Dune::XT::Common::get_matrix_rows(serializable_matrix);
#ifdef DUNE_XT_LA_DISABLE_ALL_CHECKS
  assert(Dune::XT::Common::get_matrix_cols(serializable_matrix) == size);
#else
  if (Dune::XT::Common::get_matrix_cols(serializable_matrix) != size)
    DUNE_THROW(Exceptions::eigen_solver_failed_bc_data_did_not_fulfill_requirements,
               "Given matrix has to be square, is " << size << "x"
                                                    << Dune::XT::Common::get_matrix_cols(serializable_matrix)
                                                    << "!");
#endif // DUNE_XT_LA_DISABLE_ALL_CHECKS
  std::vector<double> real_part_of_eigenvalues(size, 0.);
  std::vector<double> imag_part_of_eigenvalues(size, 0.);
  std::vector<double> dummy_left_eigenvalues(1, 0.);
  std::vector<double> dummy_right_eigenvalues(1, 0.);
  // lapacks favorite storage format is column-major, otherwise the matrix would be copied from row-major to col-major
  const int info = Common::Lapacke::dgeev(Common::Lapacke::col_major(),
                                          /*do_not_compute_left_egenvectors: */ 'N',
                                          /*do_not_compute_right_egenvectors: */ 'N',
                                          size,
                                          Dune::XT::Common::serialize_colwise<double>(serializable_matrix).get(),
                                          size,
                                          real_part_of_eigenvalues.data(),
                                          imag_part_of_eigenvalues.data(),
                                          dummy_left_eigenvalues.data(),
                                          size,
                                          dummy_right_eigenvalues.data(),
                                          size);
  if (info != 0)
    DUNE_THROW(Dune::XT::LA::Exceptions::eigen_solver_failed, "The lapack backend reported '" << info << "'!");
  std::vector<std::complex<double>> eigenvalues(size);
  for (size_t ii = 0; ii < size; ++ii)
    eigenvalues[ii] = {real_part_of_eigenvalues[ii], imag_part_of_eigenvalues[ii]};
  return eigenvalues;
} // ... compute_eigenvalues_of_a_real_matrix_using_lapack(...)


/**
 * \sa https://software.intel.com/en-us/mkl-developer-reference-c-geev
 * \note Most likely, you do not want to use this function directly, but
 *       compute_eigenvalues_and_right_eigenvectors_using_lapack.
 */
template <class SerializableRealMatrixType, class ComplexMatrixType>
typename std::enable_if<Common::is_matrix<SerializableRealMatrixType>::value
                            && Common::is_matrix<ComplexMatrixType>::value,
                        void>::type
compute_eigenvalues_and_right_eigenvectors_of_a_real_matrix_using_lapack(
    const SerializableRealMatrixType& serializable_matrix,
    std::vector<std::complex<double>>& eigenvalues,
    ComplexMatrixType& right_eigenvectors)
{
  if (!Common::Lapacke::available())
    DUNE_THROW(Exceptions::eigen_solver_failed_bc_it_was_not_set_up_correctly,
               "Do not call any lapack related method if Common::Lapacke::available() is false!");
  using real_type = typename Dune::XT::Common::MatrixAbstraction<SerializableRealMatrixType>::S;
  static_assert(Dune::XT::Common::is_arithmetic<real_type>::value && !Dune::XT::Common::is_complex<real_type>::value,
                "Not implemented for complex matrices (yet)!");
  using complex_type = typename Dune::XT::Common::MatrixAbstraction<ComplexMatrixType>::S;
  static_assert(Dune::XT::Common::is_complex<complex_type>::value,
                "You have to manually convert the eigenvaluematrix to something else outside!");
  static_assert(std::is_same<Dune::XT::Common::real_t<complex_type>, double>::value,
                "You have to manually convert the eigenvaluematrix to something else outside!");
  const size_t size = Dune::XT::Common::get_matrix_rows(serializable_matrix);
#ifdef DUNE_XT_LA_DISABLE_ALL_CHECKS
  assert(Dune::XT::Common::get_matrix_cols(serializable_matrix) == size);
  assert(Dune::XT::Common::get_matrix_rows(right_eigenvectors) == size);
  assert(Dune::XT::Common::get_matrix_cols(right_eigenvectors) == size);
#else
  if (Dune::XT::Common::get_matrix_cols(serializable_matrix) != size)
    DUNE_THROW(Exceptions::eigen_solver_failed_bc_data_did_not_fulfill_requirements,
               "Given matrix has to be square, is " << size << "x"
                                                    << Dune::XT::Common::get_matrix_cols(serializable_matrix)
                                                    << "!");
  if (Dune::XT::Common::get_matrix_rows(right_eigenvectors) != size)
    DUNE_THROW(Exceptions::eigen_solver_failed_bc_data_did_not_fulfill_requirements,
               "Given matrix of right eigenvectors has to be of same size as given matrix, is "
                   << Dune::XT::Common::get_matrix_rows(right_eigenvectors)
                   << "x"
                   << Dune::XT::Common::get_matrix_cols(right_eigenvectors)
                   << "!");
  if (Dune::XT::Common::get_matrix_cols(right_eigenvectors) != size)
    DUNE_THROW(Exceptions::eigen_solver_failed_bc_data_did_not_fulfill_requirements,
               "Given matrix of right eigenvectors has to be of same size as given matrix, is "
                   << Dune::XT::Common::get_matrix_rows(right_eigenvectors)
                   << "x"
                   << Dune::XT::Common::get_matrix_cols(right_eigenvectors)
                   << "!");
#endif // DUNE_XT_LA_DISABLE_ALL_CHECKS
  std::vector<double> real_part_of_eigenvalues(size, 0.);
  std::vector<double> imag_part_of_eigenvalues(size, 0.);
  std::vector<double> dummy_left_eigenvalues(1, 0.);
  std::vector<double> right_eigenvalues(size * size, 0.);
  // lapacks favorite storage format is column-major, otherwise the matrix would be copied from row-major to col-major
  const int info = Common::Lapacke::dgeev(Common::Lapacke::col_major(),
                                          /*do_not_compute_left_egenvectors: */ 'N',
                                          /*compute_right_egenvectors: */ 'V',
                                          size,
                                          Dune::XT::Common::serialize_colwise<double>(serializable_matrix).get(),
                                          size,
                                          real_part_of_eigenvalues.data(),
                                          imag_part_of_eigenvalues.data(),
                                          dummy_left_eigenvalues.data(),
                                          size,
                                          right_eigenvalues.data(),
                                          size);
  if (info != 0)
    DUNE_THROW(Dune::XT::LA::Exceptions::eigen_solver_failed, "The lapack backend reported '" << info << "'!");
  // set eigenvalues
  if (eigenvalues.size() != size)
    eigenvalues.resize(size);
  for (size_t ii = 0; ii < size; ++ii)
    eigenvalues[ii] = {real_part_of_eigenvalues[ii], imag_part_of_eigenvalues[ii]};
  // set eigenvectors, these are in column-major format, hence the jj/kk switch on assignment
  size_t jj = 0;
  while (jj < size) {
    const auto& imag_part_of_jth_eigenvalue = imag_part_of_eigenvalues[jj];
    if (!(imag_part_of_jth_eigenvalue < 0. || imag_part_of_jth_eigenvalue > 0.)) {
      // this is a real eigenvalue with corresponding real eigenvector
      for (size_t kk = 0; kk < size; ++kk) {
        const double kth_component_of_jth_eigenvector = right_eigenvalues[kk + jj * size];
        Dune::XT::Common::set_matrix_entry(
            right_eigenvectors, kk, jj, complex_type(kth_component_of_jth_eigenvector, 0.));
      }
      jj += 1;
    } else {
      // this eigenvalue and the next form a complex conjugate pair
      assert(jj + 1 < size && "This must not happen, the lapack documentation promised otherwise!");
      for (size_t kk = 0; kk < size; ++kk) {
        const double real_part_of_kth_component_of_jth_eigenvector = right_eigenvalues[kk + jj * size];
        const double imag_part_of_kth_component_of_jth_eigenvector = right_eigenvalues[kk + (jj + 1) * size];
        Dune::XT::Common::set_matrix_entry(
            right_eigenvectors,
            kk,
            jj,
            complex_type(real_part_of_kth_component_of_jth_eigenvector, imag_part_of_kth_component_of_jth_eigenvector));
        const double real_part_of_kth_component_of_jplusoneth_eigenvector = right_eigenvalues[kk + jj * size];
        const double imag_part_of_kth_component_of_jplusoneth_eigenvector =
            -1 * right_eigenvalues[kk + (jj + 1) * size];
        Dune::XT::Common::set_matrix_entry(right_eigenvectors,
                                           kk,
                                           jj + 1,
                                           complex_type(real_part_of_kth_component_of_jplusoneth_eigenvector,
                                                        imag_part_of_kth_component_of_jplusoneth_eigenvector));
      }
      jj += 2;
    }
  }
} // ... compute_eigenvalues_and_right_eigenvectors_of_a_real_matrix_using_lapack(...)


/**
 * \note Most likely, you do not want to use this function directly, but compute_eigenvalues_using_lapack or
 *       compute_eigenvalues_and_right_eigenvectors_using_lapack.
 */
template <class MatrixType>
struct lapack_helper
{
  static_assert(Common::is_matrix<MatrixType>::value, "");

  template <bool is_complex = Common::is_complex<typename Common::MatrixAbstraction<MatrixType>::S>::value,
            bool anything = true>
  struct dtype_switch;

  template <bool anything>
  struct dtype_switch<true, anything>
  {
    static inline std::vector<std::complex<double>> eigenvalues(const MatrixType& /*matrix*/)
    {
      static_assert(AlwaysFalse<MatrixType>::value,
                    "Not yet implemented for complex matrices, take a look at "
                    "https://software.intel.com/en-us/mkl-developer-reference-c-geev "
                    "and add a corresponding free function like "
                    "compute_eigenvalues_of_a_real_matrix_using_lapack(...)!");
      return std::vector<std::complex<double>>();
    }

    template <class V, class E>
    static inline void eigenvectors(const MatrixType& /*matrix*/, V& /*eigenvalues*/, E& /*eigenvectors*/)
    {
      static_assert(AlwaysFalse<MatrixType>::value,
                    "Not yet implemented for complex matrices, take a look at "
                    "https://software.intel.com/en-us/mkl-developer-reference-c-geev and add a corresponding free "
                    "function like compute_eigenvalues_and_right_eigenvectors_of_a_real_matrix_using_lapack(...)!");
    }
  };

  template <bool anything>
  struct dtype_switch<false, anything>
  {
    static inline std::vector<std::complex<double>> eigenvalues(const MatrixType& matrix)
    {
      return compute_eigenvalues_of_a_real_matrix_using_lapack(matrix);
    }

    template <class V, class E>
    static inline void eigenvectors(const MatrixType& matrix, V& eigenvalues, E& eigenvectors)
    {
      compute_eigenvalues_and_right_eigenvectors_of_a_real_matrix_using_lapack(matrix, eigenvalues, eigenvectors);
    }
  };
}; // class lapack_helper


template <class SerializableMatrixType>
typename std::enable_if<Common::is_matrix<SerializableMatrixType>::value, std::vector<std::complex<double>>>::type
compute_eigenvalues_using_lapack(const SerializableMatrixType& serializable_matrix)
{
  return lapack_helper<SerializableMatrixType>::template dtype_switch<>::eigenvalues(serializable_matrix);
}


template <class SerializableMatrixType, class ComplexMatrixType>
typename std::enable_if<Common::is_matrix<SerializableMatrixType>::value && Common::is_matrix<ComplexMatrixType>::value,
                        void>::type
compute_eigenvalues_and_right_eigenvectors_using_lapack(const SerializableMatrixType& serializable_matrix,
                                                        std::vector<std::complex<double>>& eigenvalues,
                                                        ComplexMatrixType& right_eigenvectors)
{
  lapack_helper<SerializableMatrixType>::template dtype_switch<>::eigenvectors(
      serializable_matrix, eigenvalues, right_eigenvectors);
}


} // namespace internal
} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_EIGEN_SOLVER_INTERNAL_LAPACKE_HH
