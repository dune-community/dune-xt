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

#ifndef DUNE_XT_LA_CONTAINER_COMMON_MATRIX_DENSE_HH
#define DUNE_XT_LA_CONTAINER_COMMON_MATRIX_DENSE_HH

#include <cmath>
#include <memory>
#include <mutex>
#include <numeric>
#include <vector>

#include <boost/align/aligned_allocator.hpp>

#include <dune/common/ftraits.hh>
#include <dune/common/unused.hh>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/vector.hh>

#include <dune/xt/la/container/matrix-interface.hh>
#include <dune/xt/la/container/pattern.hh>
#include <dune/xt/la/container/vector-view.hh>

namespace Dune {
namespace XT {
namespace LA {


template <class ScalarImp, Common::StorageLayout storage_layout>
class CommonDenseMatrix;

template <class ScalarImp>
class CommonSparseVector;


namespace internal {


template <class ScalarType>
struct MatrixBackendBase
{
  using VectorType = std::vector<ScalarType, boost::alignment::aligned_allocator<ScalarType, 64>>;
  MatrixBackendBase(const size_t num_rows, const size_t num_cols, const ScalarType value)
    : num_rows_(num_rows)
    , num_cols_(num_cols)
    , entries_(num_rows_ * num_cols_, value)
  {}

  void resize(const size_t num_rows, const size_t num_cols)
  {
    num_rows_ = num_rows;
    num_cols_ = num_cols;
    entries_.resize(num_rows_ * num_cols_);
  }

  size_t num_rows_;
  size_t num_cols_;
  VectorType entries_;
};

template <class ScalarType, Common::StorageLayout = Common::StorageLayout::dense_row_major>
struct CommonDenseMatrixBackend;

template <class ScalarType>
struct CommonDenseMatrixBackend<ScalarType, Common::StorageLayout::dense_row_major>
  : public MatrixBackendBase<ScalarType>
{
  using BaseType = MatrixBackendBase<ScalarType>;
  using typename BaseType::VectorType;

  CommonDenseMatrixBackend(const size_t num_rows, const size_t num_cols, const ScalarType value = ScalarType(0))
    : BaseType(num_rows, num_cols, value)
  {}

  ScalarType& get_entry_ref(const size_t rr, const size_t cc)
  {
    return entries_[rr * num_cols_ + cc];
  }

  const ScalarType& get_entry_ref(const size_t rr, const size_t cc) const
  {
    return entries_[rr * num_cols_ + cc];
  }

  VectorView<VectorType> get_view(const size_t row)
  {
    return VectorView<VectorType>(entries_, row * num_cols_, (row + 1) * num_cols_);
  }

  ConstVectorView<VectorType> get_const_view(const size_t row) const
  {
    return ConstVectorView<VectorType>(entries_, row * num_cols_, (row + 1) * num_cols_);
  }

  using BaseType::entries_;
  using BaseType::num_cols_;
};


template <class ScalarType>
struct CommonDenseMatrixBackend<ScalarType, Common::StorageLayout::dense_column_major>
  : public MatrixBackendBase<ScalarType>
{
  using BaseType = MatrixBackendBase<ScalarType>;
  using typename BaseType::VectorType;

  CommonDenseMatrixBackend(const size_t num_rows, const size_t num_cols, const ScalarType value = ScalarType(0))
    : BaseType(num_rows, num_cols, value)
  {}

  ScalarType& get_entry_ref(const size_t rr, const size_t cc)
  {
    return entries_[cc * num_rows_ + rr];
  }

  const ScalarType& get_entry_ref(const size_t rr, const size_t cc) const
  {
    return entries_[cc * num_rows_ + rr];
  }

  VectorView<VectorType> get_view(const size_t col)
  {
    return VectorView<VectorType>(entries_, col * num_rows_, (col + 1) * num_rows_);
  }

  ConstVectorView<VectorType> get_const_view(const size_t col) const
  {
    return ConstVectorView<VectorType>(entries_, col * num_rows_, (col + 1) * num_rows_);
  }

  using BaseType::entries_;
  using BaseType::num_rows_;
};


template <class ScalarImp, Common::StorageLayout storage_layout>
class CommonDenseMatrixTraits
  : public MatrixTraitsBase<ScalarImp,
                            CommonDenseMatrix<ScalarImp, storage_layout>,
                            CommonDenseMatrixBackend<ScalarImp, storage_layout>,
                            Backends::common_dense,
                            Backends::common_dense>
{};


} // namespace internal


/**
 *  \brief  A dense matrix implementation of MatrixInterface using a std::vector backend.
 */
template <class ScalarImp = double, Common::StorageLayout storage_layout = Common::StorageLayout::dense_row_major>
class CommonDenseMatrix
  : public MatrixInterface<internal::CommonDenseMatrixTraits<ScalarImp, storage_layout>, ScalarImp>
  , public ProvidesBackend<internal::CommonDenseMatrixTraits<ScalarImp, storage_layout>>
{
  using ThisType = CommonDenseMatrix;
  using InterfaceType = MatrixInterface<internal::CommonDenseMatrixTraits<ScalarImp, storage_layout>, ScalarImp>;

public:
  using typename InterfaceType::RealType;
  using typename InterfaceType::ScalarType;
  using Traits = typename InterfaceType::Traits;
  using typename ProvidesBackend<Traits>::BackendType;

private:
  using MutexesType = typename Traits::MutexesType;

public:
  explicit CommonDenseMatrix(const size_t rr = 0,
                             const size_t cc = 0,
                             const ScalarType value = ScalarType(0),
                             const size_t num_mutexes = 1)
    : backend_(std::make_unique<BackendType>(rr, cc, value))
    , mutexes_(std::make_unique<MutexesType>(num_mutexes))
  {}

  /// This constructors ignores the given pattern and initializes the matrix with 0.
  CommonDenseMatrix(const size_t rr,
                    const size_t cc,
                    const SparsityPatternDefault& /*pattern*/,
                    const size_t num_mutexes = 1)
    : backend_(std::make_unique<BackendType>(rr, cc, ScalarType(0)))
    , mutexes_(std::make_unique<MutexesType>(num_mutexes))
  {}

  CommonDenseMatrix(const ThisType& other)
    : backend_(std::make_unique<BackendType>(*other.backend_))
    , mutexes_(std::make_unique<MutexesType>(other.mutexes_->size()))
  {}

  /**
   * \note If prune == true, this implementation is not optimal!
   */
  explicit CommonDenseMatrix(const BackendType& other,
                             const bool prune = false,
                             const typename Common::FloatCmp::DefaultEpsilon<ScalarType>::Type eps =
                                 Common::FloatCmp::DefaultEpsilon<ScalarType>::value(),
                             const size_t num_mutexes = 1)
    : mutexes_(std::make_unique<MutexesType>(num_mutexes))
  {
    if (prune)
      backend_ = ThisType(other).pruned(eps).backend_;
    else
      backend_ = std::make_unique<BackendType>(other);
  }

  template <class T>
  CommonDenseMatrix(const DenseMatrix<T>& other, const size_t num_mutexes = 1)
    : backend_(std::make_unique<BackendType>(other.rows(), other.cols()))
    , mutexes_(std::make_unique<MutexesType>(num_mutexes))
  {
    for (size_t ii = 0; ii < other.rows(); ++ii)
      for (size_t jj = 0; jj < other.cols(); ++jj)
        set_entry(ii, jj, other[ii][jj]);
  } // CommonDenseMatrix(...)

  template <class T, class F>
  CommonDenseMatrix(const MatrixInterface<T, F>& other, const size_t num_mutexes = 1)
    : backend_(std::make_unique<BackendType>(other.rows(), other.cols()))
    , mutexes_(std::make_unique<MutexesType>(num_mutexes))
  {
    for (size_t ii = 0; ii < other.rows(); ++ii)
      for (size_t jj = 0; jj < other.cols(); ++jj)
        set_entry(ii, jj, other.get_entry(ii, jj));
  } // CommonDenseMatrix(...)

  /**
   *  \note Takes ownership of backend_ptr in the sense that you must not delete it afterwards!
   */
  explicit CommonDenseMatrix(BackendType* backend_ptr, const size_t num_mutexes = 1)
    : backend_(backend_ptr)
    , mutexes_(std::make_unique<MutexesType>(num_mutexes))
  {}

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      *backend_ = *other.backend_;
      mutexes_ = std::make_unique<MutexesType>(other.mutexes_->size());
    }
    return *this;
  }

  ThisType& copy_backend(const ThisType& other)
  {
    *backend_ = *other.backend_;
  }

  ThisType& operator=(const BackendType& other)
  {
    backend_ = std::make_unique<BackendType>(other);
    return *this;
  }

  /// \name Required by the ProvidesBackend interface.
  /// \{

  BackendType& backend()
  {
    return *backend_;
  }

  const BackendType& backend() const
  {
    return *backend_;
  }

  ScalarType* data()
  {
    return backend().entries_.data();
  }

  const ScalarType* data() const
  {
    return backend().entries_.data();
  }

  void resize(const size_t ii, const size_t jj)
  {
    backend_->resize(ii, jj);
  }

  /// \}
  /// \name Required by ContainerInterface.
  /// \{

  ThisType copy() const
  {
    return ThisType(*backend_);
  }

  void scal(const ScalarType& alpha)
  {
    const internal::VectorLockGuard DUNE_UNUSED(guard)(*mutexes_);
    for (auto& entry : backend_->entries_)
      entry *= alpha;
  }

  template <class OtherMatrixType>
  std::enable_if_t<Common::is_matrix<OtherMatrixType>::value, void> axpy(const ScalarType& alpha,
                                                                         const OtherMatrixType& xx)
  {
    if (!has_equal_shape(xx))
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The shape of xx (" << xx.rows() << "x" << xx.cols() << ") does not match the shape of this ("
                                     << rows() << "x" << cols() << ")!");
    const internal::VectorLockGuard DUNE_UNUSED(guard)(*mutexes_);
    axpy_impl(alpha, xx);
  } // ... axpy(...)

private:
  void axpy_impl(const ScalarType& alpha, const ThisType& xx)
  {
    for (size_t ii = 0; ii < backend_->entries_.size(); ++ii)
      backend_->entries_[ii] += alpha * xx.backend_->entries_[ii];
  }

  template <class OtherMatrixType>
  void axpy_impl(const ScalarType& alpha, const OtherMatrixType& xx)
  {
    for (size_t ii = 0; ii < rows(); ++ii)
      for (size_t jj = 0; jj < cols(); ++jj)
        backend_->get_entry_ref(ii, jj) += alpha * Common::MatrixAbstraction<OtherMatrixType>::get_entry(xx, ii, jj);
  }

public:
  template <class OtherMatrixType>
  std::enable_if_t<Common::is_matrix<OtherMatrixType>::value, bool> has_equal_shape(const OtherMatrixType& other) const
  {
    return (rows() == other.rows()) && (cols() == other.cols());
  }

  /// \}
  /// \name Required by MatrixInterface.
  /// \{

  inline size_t rows() const
  {
    return backend_->num_rows_;
  }

  inline size_t cols() const
  {
    return backend_->num_cols_;
  }

  template <class FirstVectorType, class SecondVectorType>
  inline void mv(const FirstVectorType& xx, SecondVectorType& yy) const
  {
    using V1 = typename Common::VectorAbstraction<FirstVectorType>;
    using V2 = typename Common::VectorAbstraction<SecondVectorType>;
    static_assert(V1::is_vector && V2::is_vector, "");
    assert(xx.size() == cols() && yy.size() == rows());
    if (storage_layout == Common::StorageLayout::dense_row_major && V1::is_contiguous) {
      for (size_t rr = 0; rr < rows(); ++rr)
        V2::set_entry(
            yy,
            rr,
            std::inner_product(&get_entry_ref(rr, 0.), &get_entry_ref(rr + 1, 0.), V1::data(xx), ScalarType(0.)));
    } else {
      assert(xx.size() == cols() && yy.size() == rows());
      yy *= ScalarType(0.);
      for (size_t rr = 0; rr < rows(); ++rr) {
        V2::set_entry(yy, rr, 0.);
        for (size_t cc = 0; cc < cols(); ++cc)
          V2::add_to_entry(yy, rr, get_entry(rr, cc) * V1::get_entry(xx, cc));
      }
    }
  }

  template <class FirstVectorType, class SecondVectorType>
  inline void mtv(const FirstVectorType& xx, SecondVectorType& yy) const
  {
    using V1 = typename Common::VectorAbstraction<FirstVectorType>;
    using V2 = typename Common::VectorAbstraction<SecondVectorType>;
    static_assert(V1::is_vector && V2::is_vector, "");
    assert(xx.size() == rows() && yy.size() == cols());
    yy *= ScalarType(0.);
    for (size_t cc = 0; cc < cols(); ++cc) {
      V2::set_entry(yy, cc, 0.);
      for (size_t rr = 0; rr < rows(); ++rr)
        V2::add_to_entry(yy, cc, get_entry(cc, rr) * V1::get_entry(xx, rr));
    }
  }

  void mtv(const CommonSparseVector<ScalarType>& xx, CommonSparseVector<ScalarType>& yy) const
  {
    yy.clear();
    const auto& vec_entries = xx.entries();
    const auto& vec_indices = xx.indices();
    thread_local std::vector<ScalarType> tmp_vec;
    tmp_vec.resize(cols(), 0.);
    std::fill(tmp_vec.begin(), tmp_vec.end(), 0.);
    for (size_t ii = 0; ii < vec_entries.size(); ++ii) {
      const size_t rr = vec_indices[ii];
      for (size_t cc = 0; cc < cols(); ++cc)
        tmp_vec[cc] += vec_entries[ii] * backend_->get_entry_ref(rr, cc);
    }
    for (size_t cc = 0; cc < cols(); ++cc)
      if (XT::Common::FloatCmp::ne(tmp_vec[cc], ScalarType(0.)))
        yy.set_new_entry(cc, tmp_vec[cc]);
  } // void mtv(...)

  void add_to_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    assert(ii < rows());
    assert(jj < cols());
    internal::LockGuard DUNE_UNUSED(lock)(*mutexes_, ii, rows());
    backend_->get_entry_ref(ii, jj) += value;
  } // ... add_to_entry(...)

  void set_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    assert(ii < rows());
    assert(jj < cols());
    backend_->get_entry_ref(ii, jj) = value;
  } // ... set_entry(...)

  void set_all_entries(const ScalarType& value)
  {
    std::fill(backend_->entries_.begin(), backend_->entries_.end(), value);
  } // ... set_entry(...)

  ScalarType get_entry(const size_t ii, const size_t jj) const
  {
    assert(ii < rows());
    assert(jj < cols());
    return backend_->get_entry_ref(ii, jj);
  } // ... get_entry(...)

  inline VectorView<typename BackendType::VectorType> operator[](const size_t ii)
  {
    return backend_->get_view(ii);
  }

  inline ConstVectorView<typename BackendType::VectorType> operator[](const size_t ii) const
  {
    return backend_->get_const_view(ii);
  }

  void clear_row(const size_t ii)
  {
    if (ii >= rows())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given ii (" << ii << ") is larger than the rows of this (" << rows() << ")!");
    std::fill_n(&(backend_->get_entry_ref(ii, 0)), cols(), ScalarType(0));
  } // ... clear_row(...)

  void clear_col(const size_t jj)
  {
    if (jj >= cols())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given jj (" << jj << ") is larger than the cols of this (" << cols() << ")!");
    for (size_t ii = 0; ii < rows(); ++ii)
      backend_->get_entry_ref(ii, jj) = ScalarType(0);
  } // ... clear_col(...)

  void unit_row(const size_t ii)
  {
    if (ii >= cols())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given ii (" << ii << ") is larger than the cols of this (" << cols() << ")!");
    if (ii >= rows())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given ii (" << ii << ") is larger than the rows of this (" << rows() << ")!");
    clear_row(ii);
    set_entry(ii, ii, ScalarType(1));
  } // ... unit_row(...)

  void unit_col(const size_t jj)
  {
    if (jj >= cols())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given jj (" << jj << ") is larger than the cols of this (" << cols() << ")!");
    if (jj >= rows())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given jj (" << jj << ") is larger than the rows of this (" << rows() << ")!");
    clear_col(jj);
    set_entry(jj, jj, ScalarType(1));
  } // ... unit_col(...)

  bool valid() const
  {
    for (const auto& entry : backend_->entries_)
      if (Common::isnan(entry) || Common::isinf(entry))
        return false;
    return true;
  } // ... valid(...)

  /// \}

  using InterfaceType::operator+;
  using InterfaceType::operator-;
  using InterfaceType::operator+=;
  using InterfaceType::operator-=;

  template <class OtherMatrixType>
  void rightmultiply(const OtherMatrixType& other)
  {
    using M = typename Common::MatrixAbstraction<OtherMatrixType>;
    static_assert(M::is_matrix, "");
    BackendType new_backend(rows(), M::cols(other), ScalarType(0.));
    if (M::rows(other) != cols())
      DUNE_THROW(Dune::XT::Common::Exceptions::shapes_do_not_match,
                 "For rightmultiply, the number of columns of this has to match the number of rows of other!");
    for (size_t rr = 0; rr < rows(); ++rr)
      for (size_t cc = 0; cc < cols(); ++cc)
        for (size_t kk = 0; kk < cols(); ++kk)
          new_backend->get_entry_ref(rr, cc) += get_entry(rr, kk) * M::get_entry(other, kk, cc);
    *backend_ = new_backend;
  }

  virtual ThisType pruned(const typename Common::FloatCmp::DefaultEpsilon<ScalarType>::Type eps =
                              Common::FloatCmp::DefaultEpsilon<ScalarType>::value()) const override
  {
    auto ret = this->copy();
    for (auto& entry : ret.backend_->entries_)
      if (XT::Common::FloatCmp::eq<Common::FloatCmp::Style::absolute>(ScalarType(0.), entry, 0., eps))
        entry = ScalarType(0);
    return ret;
  } // ... pruned(...)

  ScalarType& get_entry_ref(const size_t rr, const size_t cc)
  {
    return backend_->get_entry_ref(rr, cc);
  }

  const ScalarType& get_entry_ref(const size_t rr, const size_t cc) const
  {
    return backend_->get_entry_ref(rr, cc);
  }

  // get pointer to entry
  ScalarType* get_entry_ptr(const size_t rr, const size_t cc)
  {
    return &(backend_->get_entry_ref(rr, cc));
  }

private:
  std::unique_ptr<BackendType> backend_;
  std::unique_ptr<MutexesType> mutexes_;
}; // class CommonDenseMatrix


} // namespace LA
namespace Common {


template <class T, Common::StorageLayout layout>
struct MatrixAbstraction<LA::CommonDenseMatrix<T, layout>>
  : public LA::internal::MatrixAbstractionBase<LA::CommonDenseMatrix<T, layout>>
{
  using BaseType = LA::internal::MatrixAbstractionBase<LA::CommonDenseMatrix<T, layout>>;

  static const constexpr Common::StorageLayout storage_layout = layout;

  template <size_t rows = BaseType::static_rows, size_t cols = BaseType::static_cols, class FieldType = T>
  using MatrixTypeTemplate = LA::CommonDenseMatrix<FieldType, layout>;

  static inline T* data(LA::CommonDenseMatrix<T, layout>& mat)
  {
    return mat.data();
  }

  static inline const T* data(const LA::CommonDenseMatrix<T, layout>& mat)
  {
    return mat.data();
  }
};


} // namespace Common
} // namespace XT
} // namespace Dune


// begin: this is what we need for the lib
#if DUNE_XT_WITH_PYTHON_BINDINGS


extern template class Dune::XT::LA::CommonDenseMatrix<double>;


#endif // DUNE_XT_WITH_PYTHON_BINDINGS
// end: this is what we need for the lib


#endif // DUNE_XT_LA_CONTAINER_COMMON_MATRIX_DENSE_HH
