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

#include <dune/xt/la/container/common/matrix/dense.hh>
#include <dune/xt/la/exceptions.hh>
#include <dune/xt/la/type_traits.hh>

namespace Dune {
namespace XT {
namespace LA {
namespace internal {


template <class MatrixType>
struct is_contiguous_and_mutable
{
  static constexpr bool value =
      (Common::MatrixAbstraction<std::decay_t<MatrixType>>::storage_layout == Common::StorageLayout::dense_row_major
       || Common::MatrixAbstraction<std::decay_t<MatrixType>>::storage_layout
              == Common::StorageLayout::dense_column_major)
      && !std::is_const<std::remove_reference_t<MatrixType>>::value;
};


template <class MatrixType, bool contiguous_and_mutable>
class MatrixDataProvider;

template <class MatrixType>
class MatrixDataProvider<MatrixType, true>
{
public:
  MatrixDataProvider(MatrixType& matrix)
    : matrix_(matrix)
  {}

  double* data()
  {
    return Common::MatrixAbstraction<MatrixType>::data(matrix_);
  }

  static const Common::StorageLayout storage_layout = Common::MatrixAbstraction<MatrixType>::storage_layout;

private:
  MatrixType& matrix_;
};

template <class MatrixType>
class MatrixDataProvider<MatrixType, false>
{
public:
  // lapacks favorite storage format is column-major, otherwise the matrix would be copied from row-major to col-major
  MatrixDataProvider(const MatrixType& matrix)
    : serialized_matrix_(Dune::XT::Common::serialize_colwise<double>(matrix))
  {}

  double* data()
  {
    return serialized_matrix_.get();
  }

  static const Common::StorageLayout storage_layout = Common::StorageLayout::dense_column_major;

private:
  std::unique_ptr<double[]> serialized_matrix_;
};


/**
 * \sa https://software.intel.com/en-us/mkl-developer-reference-c-geev
 * \note Most likely, you do not want to use this function directly, but compute_eigenvalues_using_lapack.
 */
template <class RealMatrixType>
typename std::enable_if<Common::is_matrix<std::decay_t<RealMatrixType>>::value, std::vector<std::complex<double>>>::type
compute_eigenvalues_of_a_real_matrix_using_lapack(RealMatrixType&& matrix)
{
  MatrixDataProvider<std::decay_t<RealMatrixType>, is_contiguous_and_mutable<RealMatrixType>::value>
      matrix_data_provider(matrix);
  return compute_eigenvalues_of_a_real_matrix_using_lapack_impl(matrix, matrix_data_provider);
}


template <class RealMatrixType, bool contiguous_and_mutable>
typename std::enable_if<Common::is_matrix<RealMatrixType>::value, std::vector<std::complex<double>>>::type
compute_eigenvalues_of_a_real_matrix_using_lapack_impl(
    const RealMatrixType& matrix, MatrixDataProvider<RealMatrixType, contiguous_and_mutable>& matrix_data_provider)
{
  if (!Common::Lapacke::available())
    DUNE_THROW(Exceptions::eigen_solver_failed_bc_it_was_not_set_up_correctly,
               "Do not call any lapack related method if Common::Lapacke::available() is false!");
  using real_type = typename Dune::XT::Common::MatrixAbstraction<RealMatrixType>::S;
  static_assert(Dune::XT::Common::is_arithmetic<real_type>::value && !Dune::XT::Common::is_complex<real_type>::value,
                "Not implemented for complex matrices (yet)!");
  const size_t size = Dune::XT::Common::get_matrix_rows(matrix);
#ifdef DUNE_XT_LA_DISABLE_ALL_CHECKS
  assert(Dune::XT::Common::get_matrix_cols(matrix) == size);
#else
  if (Dune::XT::Common::get_matrix_cols(matrix) != size)
    DUNE_THROW(Exceptions::eigen_solver_failed_bc_data_did_not_fulfill_requirements,
               "Given matrix has to be square, is " << size << "x" << Dune::XT::Common::get_matrix_cols(matrix) << "!");
#endif // DUNE_XT_LA_DISABLE_ALL_CHECKS
  thread_local std::vector<double> real_part_of_eigenvalues(size, 0.);
  thread_local std::vector<double> imag_part_of_eigenvalues(size, 0.);
  assert(size < std::numeric_limits<int>::max());
  int storage_layout = MatrixDataProvider<RealMatrixType, contiguous_and_mutable>::storage_layout
                               == Common::StorageLayout::dense_row_major
                           ? Common::Lapacke::row_major()
                           : Common::Lapacke::col_major();
  thread_local std::vector<double> work(1);
  thread_local size_t last_size = -1;
  if (size != last_size) {
    real_part_of_eigenvalues.resize(size);
    imag_part_of_eigenvalues.resize(size);
    // get optimal working size in work[0] (requested by lwork = -1)
    const int info = Common::Lapacke::dgeev_work(storage_layout,
                                                 /*do_not_compute_left_eigenvectors: */ 'N',
                                                 /*do_not_compute_right_eigenvectors: */ 'N',
                                                 static_cast<int>(size),
                                                 matrix_data_provider.data(),
                                                 static_cast<int>(size),
                                                 real_part_of_eigenvalues.data(),
                                                 imag_part_of_eigenvalues.data(),
                                                 nullptr,
                                                 static_cast<int>(size),
                                                 nullptr,
                                                 static_cast<int>(size),
                                                 work.data(),
                                                 -1);
    if (info != 0)
      DUNE_THROW(Dune::XT::LA::Exceptions::eigen_solver_failed, "The lapack backend reported '" << info << "'!");
    work.resize(work[0]);
    last_size = size;
  }
  const int info = Common::Lapacke::dgeev_work(storage_layout,
                                               /*do_not_compute_left_eigenvectors: */ 'N',
                                               /*do_not_compute_right_eigenvectors: */ 'N',
                                               static_cast<int>(size),
                                               matrix_data_provider.data(),
                                               static_cast<int>(size),
                                               real_part_of_eigenvalues.data(),
                                               imag_part_of_eigenvalues.data(),
                                               nullptr,
                                               static_cast<int>(size),
                                               nullptr,
                                               static_cast<int>(size),
                                               work.data(),
                                               static_cast<int>(work.size()));
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
template <class RealMatrixType, class ComplexMatrixType>
typename std::enable_if<Common::is_matrix<std::decay_t<RealMatrixType>>::value
                            && Common::is_matrix<ComplexMatrixType>::value,
                        void>::type
compute_eigenvalues_and_right_eigenvectors_of_a_real_matrix_using_lapack(RealMatrixType&& matrix,
                                                                         std::vector<std::complex<double>>& eigenvalues,
                                                                         ComplexMatrixType& right_eigenvectors)
{
  MatrixDataProvider<std::decay_t<RealMatrixType>, is_contiguous_and_mutable<RealMatrixType>::value>
      matrix_data_provider(matrix);
  compute_eigenvalues_and_right_eigenvectors_of_a_real_matrix_using_lapack_impl(
      matrix, matrix_data_provider, eigenvalues, right_eigenvectors);
}


template <class RealMatrixType, class ComplexMatrixType, bool contiguous_and_mutable>
typename std::enable_if<Common::is_matrix<RealMatrixType>::value && Common::is_matrix<ComplexMatrixType>::value,
                        void>::type
compute_eigenvalues_and_right_eigenvectors_of_a_real_matrix_using_lapack_impl(
    const RealMatrixType& matrix,
    MatrixDataProvider<RealMatrixType, contiguous_and_mutable>& matrix_data_provider,
    std::vector<std::complex<double>>& eigenvalues,
    ComplexMatrixType& right_eigenvectors)
{
  if (!Common::Lapacke::available())
    DUNE_THROW(Exceptions::eigen_solver_failed_bc_it_was_not_set_up_correctly,
               "Do not call any lapack related method if Common::Lapacke::available() is false!");
  using real_type = typename Dune::XT::Common::MatrixAbstraction<RealMatrixType>::S;
  static_assert(Dune::XT::Common::is_arithmetic<real_type>::value && !Dune::XT::Common::is_complex<real_type>::value,
                "Not implemented for complex matrices (yet)!");
  using complex_type = typename Dune::XT::Common::MatrixAbstraction<ComplexMatrixType>::S;
  static_assert(Dune::XT::Common::is_complex<complex_type>::value,
                "You have to manually convert the eigenvaluematrix to something else outside!");
  static_assert(std::is_same<Dune::XT::Common::real_t<complex_type>, double>::value,
                "You have to manually convert the eigenvaluematrix to something else outside!");
  const size_t size = Dune::XT::Common::get_matrix_rows(matrix);
#ifdef DUNE_XT_LA_DISABLE_ALL_CHECKS
  assert(Dune::XT::Common::get_matrix_cols(matrix) == size);
  assert(Dune::XT::Common::get_matrix_rows(right_eigenvectors) == size);
  assert(Dune::XT::Common::get_matrix_cols(right_eigenvectors) == size);
#else
  if (Dune::XT::Common::get_matrix_cols(matrix) != size)
    DUNE_THROW(Exceptions::eigen_solver_failed_bc_data_did_not_fulfill_requirements,
               "Given matrix has to be square, is " << size << "x" << Dune::XT::Common::get_matrix_cols(matrix) << "!");
  if (Dune::XT::Common::get_matrix_rows(right_eigenvectors) != size)
    DUNE_THROW(Exceptions::eigen_solver_failed_bc_data_did_not_fulfill_requirements,
               "Given matrix of right eigenvectors has to be of same size as given matrix, is "
                   << Dune::XT::Common::get_matrix_rows(right_eigenvectors) << "x"
                   << Dune::XT::Common::get_matrix_cols(right_eigenvectors) << "!");
  if (Dune::XT::Common::get_matrix_cols(right_eigenvectors) != size)
    DUNE_THROW(Exceptions::eigen_solver_failed_bc_data_did_not_fulfill_requirements,
               "Given matrix of right eigenvectors has to be of same size as given matrix, is "
                   << Dune::XT::Common::get_matrix_rows(right_eigenvectors) << "x"
                   << Dune::XT::Common::get_matrix_cols(right_eigenvectors) << "!");
#endif // DUNE_XT_LA_DISABLE_ALL_CHECKS
  int storage_layout = MatrixDataProvider<RealMatrixType, contiguous_and_mutable>::storage_layout
                               == Common::StorageLayout::dense_row_major
                           ? Common::Lapacke::row_major()
                           : Common::Lapacke::col_major();
  thread_local std::vector<double> real_part_of_eigenvalues(size);
  thread_local std::vector<double> imag_part_of_eigenvalues(size);
  thread_local CommonDenseMatrix<double, MatrixDataProvider<RealMatrixType, contiguous_and_mutable>::storage_layout>
      right_eigenvectors_matrix(size, size, 0.);
  assert(size < std::numeric_limits<int>::max());
  thread_local std::vector<double> work(1);
  thread_local size_t last_size = -1;
  if (size != last_size) {
    real_part_of_eigenvalues.resize(size);
    imag_part_of_eigenvalues.resize(size);
    right_eigenvectors_matrix.backend().resize(size, size);
    // get optimal working size in work[0] (requested by lwork = -1)
    int info = Common::Lapacke::dgeev_work(storage_layout,
                                           /*do_not_compute_left_eigenvectors: */ 'N',
                                           /*compute_right_eigenvectors: */ 'V',
                                           static_cast<int>(size),
                                           matrix_data_provider.data(),
                                           static_cast<int>(size),
                                           real_part_of_eigenvalues.data(),
                                           imag_part_of_eigenvalues.data(),
                                           nullptr,
                                           static_cast<int>(size),
                                           right_eigenvectors_matrix.data(),
                                           static_cast<int>(size),
                                           work.data(),
                                           -1);
    if (info != 0)
      DUNE_THROW(Dune::XT::LA::Exceptions::eigen_solver_failed, "The lapack backend reported '" << info << "'!");
    work.resize(work[0]);
    last_size = size;
  }
  // do the actual calculation
  int info = Common::Lapacke::dgeev_work(storage_layout,
                                         /*do_not_compute_left_eigenvectors: */ 'N',
                                         /*compute_right_eigenvectors: */ 'V',
                                         static_cast<int>(size),
                                         matrix_data_provider.data(),
                                         static_cast<int>(size),
                                         real_part_of_eigenvalues.data(),
                                         imag_part_of_eigenvalues.data(),
                                         nullptr,
                                         static_cast<int>(size),
                                         right_eigenvectors_matrix.data(),
                                         static_cast<int>(size),
                                         work.data(),
                                         static_cast<int>(work.size()));
  if (info != 0)
    DUNE_THROW(Dune::XT::LA::Exceptions::eigen_solver_failed, "The lapack backend reported '" << info << "'!");
  // set eigenvalues
  if (eigenvalues.size() != size)
    eigenvalues.resize(size);
  for (size_t ii = 0; ii < size; ++ii)
    eigenvalues[ii] = {real_part_of_eigenvalues[ii], imag_part_of_eigenvalues[ii]};
  size_t jj = 0;
  while (jj < size) {
    const auto& imag_part_of_jth_eigenvalue = imag_part_of_eigenvalues[jj];
    if (imag_part_of_jth_eigenvalue == 0.) {
      // this is a real eigenvalue with corresponding real eigenvector
      for (size_t kk = 0; kk < size; ++kk) {
        Dune::XT::Common::set_matrix_entry(
            right_eigenvectors, kk, jj, complex_type(right_eigenvectors_matrix.get_entry(kk, jj), 0.));
      }
      jj += 1;
    } else {
      // this eigenvalue and the next form a complex conjugate pair
      assert(jj + 1 < size && "This must not happen, the lapack documentation promised otherwise!");
      for (size_t kk = 0; kk < size; ++kk) {
        const double real_part_of_kth_component_of_jth_eigenvector = right_eigenvectors_matrix.get_entry(kk, jj);
        const double imag_part_of_kth_component_of_jth_eigenvector = right_eigenvectors_matrix.get_entry(kk, jj + 1);
        Dune::XT::Common::set_matrix_entry(
            right_eigenvectors,
            kk,
            jj,
            complex_type(real_part_of_kth_component_of_jth_eigenvector, imag_part_of_kth_component_of_jth_eigenvector));
        const double real_part_of_kth_component_of_jplusoneth_eigenvector =
            real_part_of_kth_component_of_jth_eigenvector;
        const double imag_part_of_kth_component_of_jplusoneth_eigenvector =
            -1 * imag_part_of_kth_component_of_jth_eigenvector;
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
    template <class MatrixImp>
    static inline std::vector<std::complex<double>> eigenvalues(MatrixImp&& /*matrix*/)
    {
      static_assert(AlwaysFalse<MatrixImp>::value,
                    "Not yet implemented for complex matrices, take a look at "
                    "https://software.intel.com/en-us/mkl-developer-reference-c-geev "
                    "and add a corresponding free function like "
                    "compute_eigenvalues_of_a_real_matrix_using_lapack(...)!");
      return std::vector<std::complex<double>>();
    }

    template <class V, class E, class MatrixImp>
    static inline void eigenvectors(MatrixImp&& /*matrix*/, V& /*eigenvalues*/, E& /*eigenvectors*/)
    {
      static_assert(AlwaysFalse<MatrixImp>::value,
                    "Not yet implemented for complex matrices, take a look at "
                    "https://software.intel.com/en-us/mkl-developer-reference-c-geev and add a corresponding free "
                    "function like compute_eigenvalues_and_right_eigenvectors_of_a_real_matrix_using_lapack(...)!");
    }
  };

  template <bool anything>
  struct dtype_switch<false, anything>
  {
    template <class MatrixImp>
    static inline std::vector<std::complex<double>> eigenvalues(MatrixImp&& matrix)
    {
      return compute_eigenvalues_of_a_real_matrix_using_lapack(std::forward<MatrixImp>(matrix));
    }

    template <class V, class E, class MatrixImp>
    static inline void eigenvectors(MatrixImp&& matrix, V& eigenvalues, E& eigenvectors)
    {
      compute_eigenvalues_and_right_eigenvectors_of_a_real_matrix_using_lapack(
          std::forward<MatrixImp>(matrix), eigenvalues, eigenvectors);
    }
  };
}; // class lapack_helper


template <class MatrixType>
typename std::enable_if<Common::is_matrix<std::decay_t<MatrixType>>::value, std::vector<std::complex<double>>>::type
compute_eigenvalues_using_lapack(MatrixType&& matrix)
{
  return lapack_helper<std::decay_t<MatrixType>>::template dtype_switch<>::eigenvalues(
      std::forward<MatrixType>(matrix));
}


template <class MatrixType, class ComplexMatrixType>
typename std::enable_if<Common::is_matrix<std::decay_t<MatrixType>>::value
                            && Common::is_matrix<ComplexMatrixType>::value,
                        void>::type
compute_eigenvalues_and_right_eigenvectors_using_lapack(MatrixType&& matrix,
                                                        std::vector<std::complex<double>>& eigenvalues,
                                                        ComplexMatrixType& right_eigenvectors)
{
  lapack_helper<std::decay_t<MatrixType>>::template dtype_switch<>::eigenvectors(
      std::forward<MatrixType>(matrix), eigenvalues, right_eigenvectors);
}


} // namespace internal
} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_EIGEN_SOLVER_INTERNAL_LAPACKE_HH
