// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2017 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Barbara Verfürth (2015)
//   Felix Schindler  (2013 - 2017)
//   Rene Milk        (2014 - 2016)
//   Tobias Leibner   (2014, 2016)

#ifndef DUNE_XT_LA_CONTAINER_COMMON_HH
#define DUNE_XT_LA_CONTAINER_COMMON_HH

#include <cmath>
#include <initializer_list>
#include <memory>
#include <type_traits>
#include <vector>
#include <complex>
#include <mutex>

#include <dune/common/dynvector.hh>
#include <dune/common/dynmatrix.hh>
#include <dune/common/densematrix.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/ftraits.hh>
#include <dune/common/unused.hh>

#include <dune/xt/common/exceptions.hh>

#include "interfaces.hh"
#include "pattern.hh"

namespace Dune {
namespace XT {
namespace LA {

// forwards
template <class ScalarImp>
class CommonDenseVector;

template <class ScalarImp>
class CommonDenseMatrix;

template <class ScalarImp>
class CommonSparseMatrix;

namespace internal {

/// Traits for CommonDenseVector
template <class ScalarImp = double>
class CommonDenseVectorTraits
{
public:
  typedef typename Dune::FieldTraits<ScalarImp>::field_type ScalarType;
  typedef typename Dune::FieldTraits<ScalarImp>::real_type RealType;
  typedef ScalarType DataType;
  typedef CommonDenseVector<ScalarType> derived_type;
  typedef Dune::DynamicVector<ScalarType> BackendType;
  static const Backends backend_type = Backends::common_dense;
  static const Backends dense_matrix_type = Backends::common_dense;
  static const Backends sparse_matrix_type = Backends::common_sparse;
};

template <class ScalarImp = double>
class CommonDenseMatrixTraits
{
public:
  typedef typename Dune::FieldTraits<ScalarImp>::field_type ScalarType;
  typedef typename Dune::FieldTraits<ScalarImp>::real_type RealType;
  typedef CommonDenseMatrix<ScalarType> derived_type;
  typedef Dune::DynamicMatrix<ScalarType> BackendType;
  static const Backends backend_type = Backends::common_dense;
  static const Backends vector_type = Backends::common_dense;
};

template <class ScalarImp = double>
class CommonSparseMatrixTraits
{
public:
  typedef typename Dune::FieldTraits<ScalarImp>::field_type ScalarType;
  typedef typename Dune::FieldTraits<ScalarImp>::real_type RealType;
  typedef CommonSparseMatrix<ScalarType> derived_type;
  typedef std::vector<ScalarType> EntriesVectorType;
  static const Backends backend_type = Backends::common_dense;
  static const constexpr Backends vector_type = Backends::common_dense;
};

} // namespace internal

/**
 *  \brief A dense vector implementation of VectorInterface using the Dune::DynamicVector.
 */
template <class ScalarImp = double>
class CommonDenseVector : public VectorInterface<internal::CommonDenseVectorTraits<ScalarImp>, ScalarImp>,
                          public ProvidesBackend<internal::CommonDenseVectorTraits<ScalarImp>>,
                          public ProvidesDataAccess<internal::CommonDenseVectorTraits<ScalarImp>>
{
  typedef CommonDenseVector<ScalarImp> ThisType;
  typedef VectorInterface<internal::CommonDenseVectorTraits<ScalarImp>, ScalarImp> VectorInterfaceType;

public:
  typedef internal::CommonDenseVectorTraits<ScalarImp> Traits;
  typedef typename Traits::ScalarType ScalarType;
  typedef typename Traits::RealType RealType;
  typedef typename Traits::DataType DataType;
  typedef typename Traits::BackendType BackendType;

  explicit CommonDenseVector(const size_t ss = 0, const ScalarType value = ScalarType(0))
    : backend_(new BackendType(ss, value))
  {
  }

  explicit CommonDenseVector(const std::vector<ScalarType>& other)
    : backend_(new BackendType(other.size()))
  {
    for (size_t ii = 0; ii < other.size(); ++ii)
      backend_->operator[](ii) = other[ii];
  }

  explicit CommonDenseVector(const std::initializer_list<ScalarType>& other)
    : backend_(new BackendType(other.size()))
  {
    size_t ii = 0;
    for (auto element : other) {
      backend_->operator[](ii) = element;
      ++ii;
    }
  } // CommonDenseVector(...)

  CommonDenseVector(const ThisType& other)
    : backend_(other.backend_)
  {
  }

  explicit CommonDenseVector(const BackendType& other,
                             const bool /*prune*/ = false,
                             const ScalarType /*eps*/ = Common::FloatCmp::DefaultEpsilon<ScalarType>::value())
    : backend_(new BackendType(other))
  {
  }

  /**
   *  \note Takes ownership of backend_ptr in the sense that you must not delete it afterwards!
   */
  explicit CommonDenseVector(BackendType* backend_ptr)
    : backend_(backend_ptr)
  {
  }

  explicit CommonDenseVector(std::shared_ptr<BackendType> backend_ptr)
    : backend_(backend_ptr)
  {
  }

  ThisType& operator=(const ThisType& other)
  {
    backend_ = other.backend_;
    return *this;
  }

  ThisType& operator=(const ScalarType& value)
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    ensure_uniqueness();
    for (auto& element : *this)
      element = value;
    return *this;
  } // ... operator=(...)

  /**
   *  \note Does a deep copy.
   */
  ThisType& operator=(const BackendType& other)
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    backend_ = std::make_shared<BackendType>(other);
    return *this;
  }

  /// \name Required by the ProvidesBackend interface.
  /// \{

  BackendType& backend()
  {
    ensure_uniqueness();
    return *backend_;
  }

  const BackendType& backend() const
  {
    return *backend_;
  }

  /// \}
  /// \name Required by ProvidesDataAccess.
  /// \{

  /** \attention This makes only sense for scalar data types, not for complex! **/
  ScalarType* data()
  {
    return &(backend()[0]);
  }

  size_t data_size() const
  {
    return size();
  }

  /// \}
  /// \name Required by ContainerInterface.
  /// \{

  ThisType copy() const
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    return ThisType(*backend_);
  }

  void scal(const ScalarType& alpha)
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    auto& backend_ref = backend();
    backend_ref *= alpha;
  }

  void axpy(const ScalarType& alpha, const ThisType& xx)
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    auto& backend_ref = backend();
    if (xx.size() != size())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of x (" << xx.size() << ") does not match the size of this (" << size() << ")!");
    for (size_t ii = 0; ii < backend_ref.size(); ++ii)
      backend_ref[ii] += alpha * xx.backend()[ii];
  } // ... axpy(...)

  bool has_equal_shape(const ThisType& other) const
  {
    return size() == other.size();
  }

  /// \}
  /// \name Required by VectorInterface.
  /// \{

  inline size_t size() const
  {
    return backend_->size();
  }

  void add_to_entry(const size_t ii, const ScalarType& value)
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    auto& backend_ref = backend();
    assert(ii < size());
    backend_ref[ii] += value;
  }

  void set_entry(const size_t ii, const ScalarType& value)
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    auto& backend_ref = backend();
    assert(ii < size());
    backend_ref[ii] = value;
  }

  ScalarType get_entry(const size_t ii) const
  {
    assert(ii < size());
    return backend()[ii];
  }

protected:
  inline ScalarType& get_entry_ref(const size_t ii)
  {
    return backend_->operator[](ii);
  }

  inline const ScalarType& get_entry_ref(const size_t ii) const
  {
    return backend_->operator[](ii);
  }

public:
  /// \}
  /// \name These methods override default implementations from VectorInterface.
  /// \{

  virtual ScalarType dot(const ThisType& other) const override final
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    if (other.size() != size())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    return backend() * other.backend();
  } // ... dot(...)

  virtual RealType l1_norm() const override final
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    return backend().one_norm();
  }

  virtual RealType l2_norm() const override final
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    return backend().two_norm();
  }

  virtual RealType sup_norm() const override final
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    return backend().infinity_norm();
  }

  virtual void iadd(const ThisType& other) override final
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    auto& backend_ref = backend();
    if (other.size() != size())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    backend_ref += other.backend();
  } // ... iadd(...)

  virtual void isub(const ThisType& other) override final
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    auto& backend_ref = backend();
    if (other.size() != size())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    backend_ref -= other.backend();
  } // ... isub(...)

protected:
  /**
   * \see ContainerInterface
   */
  inline void ensure_uniqueness()
  {
    if (!backend_.unique())
      backend_ = std::make_shared<BackendType>(*backend_);
  } // ... ensure_uniqueness(...)

private:
  friend class VectorInterface<internal::CommonDenseVectorTraits<ScalarType>, ScalarType>;
  friend class CommonDenseMatrix<ScalarType>;

  std::shared_ptr<BackendType> backend_;
  mutable std::mutex mutex_;
}; // class CommonDenseVector


/**
 *  \brief  A dense matrix implementation of MatrixInterface using the Dune::DynamicMatrix.
 */
template <class ScalarImp = double>
class CommonDenseMatrix : public MatrixInterface<internal::CommonDenseMatrixTraits<ScalarImp>, ScalarImp>,
                          public ProvidesBackend<internal::CommonDenseMatrixTraits<ScalarImp>>
{
  typedef CommonDenseMatrix<ScalarImp> ThisType;
  typedef MatrixInterface<internal::CommonDenseMatrixTraits<ScalarImp>, ScalarImp> MatrixInterfaceType;

public:
  typedef internal::CommonDenseMatrixTraits<ScalarImp> Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::ScalarType ScalarType;
  typedef typename Traits::RealType RealType;

  explicit CommonDenseMatrix(const size_t rr = 0, const size_t cc = 0, const ScalarType value = ScalarType(0))
    : backend_(new BackendType(rr, cc, value))
  {
  }

  /// This constructors ignores the given pattern and initializes the matrix with 0.
  CommonDenseMatrix(const size_t rr, const size_t cc, const SparsityPatternDefault& /*pattern*/)
    : backend_(new BackendType(rr, cc, ScalarType(0)))
  {
  }

  CommonDenseMatrix(const ThisType& other)
    : backend_(other.backend_)
  {
  }

  /**
   * \note If prune == true, this implementation is not optimal!
   */
  explicit CommonDenseMatrix(const BackendType& other,
                             const bool prune = false,
                             const typename Common::FloatCmp::DefaultEpsilon<ScalarType>::Type eps =
                                 Common::FloatCmp::DefaultEpsilon<ScalarType>::value())
  {
    if (prune)
      backend_ = ThisType(other).pruned(eps).backend_;
    else
      backend_ = std::make_shared<BackendType>(other);
  }

  template <class T>
  CommonDenseMatrix(const DenseMatrix<T>& other)
    : backend_(new BackendType(other.rows(), other.cols()))
  {
    for (size_t ii = 0; ii < other.rows(); ++ii)
      for (size_t jj = 0; jj < other.cols(); ++jj)
        set_entry(ii, jj, other[ii][jj]);
  } // CommonDenseMatrix(...)

  /**
   *  \note Takes ownership of backend_ptr in the sense that you must not delete it afterwards!
   */
  explicit CommonDenseMatrix(BackendType* backend_ptr)
    : backend_(backend_ptr)
  {
  }

  explicit CommonDenseMatrix(std::shared_ptr<BackendType> backend_ptr)
    : backend_(backend_ptr)
  {
  }

  ThisType& operator=(const ThisType& other)
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    backend_ = other.backend_;
    return *this;
  }

  /**
   *  \note Does a deep copy.
   */
  ThisType& operator=(const BackendType& other)
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    backend_ = std::make_shared<BackendType>(other);
    return *this;
  }

  /// \name Required by the ProvidesBackend interface.
  /// \{

  BackendType& backend()
  {
    ensure_uniqueness();
    return *backend_;
  }

  const BackendType& backend() const
  {
    return *backend_;
  }

  /// \}
  /// \name Required by ContainerInterface.
  /// \{

  ThisType copy() const
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    return ThisType(*backend_);
  }

  void scal(const ScalarType& alpha)
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    auto& backend_ref = backend();
    backend_ref *= alpha;
  }

  void axpy(const ScalarType& alpha, const ThisType& xx)
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    auto& backend_ref = backend();
    if (!has_equal_shape(xx))
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The shape of xx (" << xx.rows() << "x" << xx.cols() << ") does not match the shape of this ("
                                     << rows()
                                     << "x"
                                     << cols()
                                     << ")!");
    backend_ref.axpy(alpha, xx.backend());
  } // ... axpy(...)

  bool has_equal_shape(const ThisType& other) const
  {
    return (rows() == other.rows()) && (cols() == other.cols());
  }

  /// \}
  /// \name Required by MatrixInterface.
  /// \{

  inline size_t rows() const
  {
    return backend().rows();
  }

  inline size_t cols() const
  {
    return backend().cols();
  }

  inline void mv(const VectorInterface<internal::CommonDenseVectorTraits<ScalarType>, ScalarType>& xx,
                 VectorInterface<internal::CommonDenseVectorTraits<ScalarType>, ScalarType>& yy) const
  {
    mv(xx.as_imp(), yy.as_imp());
  }

  inline void mv(const CommonDenseVector<ScalarType>& xx, CommonDenseVector<ScalarType>& yy) const
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    backend().mv(xx.backend(), yy.backend());
  }

  void add_to_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    auto& backend_ref = backend();
    assert(ii < rows());
    assert(jj < cols());
    backend_ref[ii][jj] += value;
  } // ... add_to_entry(...)

  void set_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    auto& backend_ref = backend();
    assert(ii < rows());
    assert(jj < cols());
    backend_ref[ii][jj] = value;
  } // ... set_entry(...)

  ScalarType get_entry(const size_t ii, const size_t jj) const
  {
    assert(ii < rows());
    assert(jj < cols());
    return backend()[ii][jj];
  } // ... get_entry(...)

  ScalarType& get_entry_ref(const size_t ii, const size_t jj)
  {
    assert(ii < rows());
    assert(jj < cols());
    return backend()[ii][jj];
  } // ... get_entry(...)

  void clear_row(const size_t ii)
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    auto& backend_ref = backend();
    if (ii >= rows())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given ii (" << ii << ") is larger than the rows of this (" << rows() << ")!");
    backend_ref[ii] *= ScalarType(0);
  } // ... clear_row(...)

  void clear_col(const size_t jj)
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    auto& backend_ref = backend();
    if (jj >= cols())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given jj (" << jj << ") is larger than the cols of this (" << cols() << ")!");
    for (size_t ii = 0; ii < rows(); ++ii)
      backend_ref[ii][jj] = ScalarType(0);
  } // ... clear_col(...)

  void unit_row(const size_t ii)
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    auto& backend_ref = backend();
    if (ii >= cols())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given ii (" << ii << ") is larger than the cols of this (" << cols() << ")!");
    if (ii >= rows())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given ii (" << ii << ") is larger than the rows of this (" << rows() << ")!");
    auto& row = backend_ref[ii];
    for (size_t jj = 0; jj < cols(); ++jj)
      row[jj] = ScalarType(0);
    row[ii] = ScalarType(1);
  } // ... unit_row(...)

  void unit_col(const size_t jj)
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    auto& backend_ref = backend();
    if (jj >= cols())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given jj (" << jj << ") is larger than the cols of this (" << cols() << ")!");
    if (jj >= rows())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given jj (" << jj << ") is larger than the rows of this (" << rows() << ")!");
    for (size_t ii = 0; ii < rows(); ++ii)
      backend_ref[ii][jj] = ScalarType(0);
    backend_ref[jj][jj] = ScalarType(1);
  } // ... unit_col(...)

  bool valid() const
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    for (size_t ii = 0; ii < rows(); ++ii) {
      const auto& row_vec = backend()[ii];
      for (size_t jj = 0; jj < cols(); ++jj) {
        const auto& entry = row_vec[jj];
        if (Common::isnan(entry) || Common::isinf(entry))
          return false;
      }
    }
    return true;
  } // ... valid(...)

  /// \}

protected:
  /**
   * \see ContainerInterface
   */
  inline void ensure_uniqueness()
  {
    if (!backend_.unique())
      backend_ = std::make_shared<BackendType>(*backend_);
  } // ... ensure_uniqueness(...)

private:
  std::shared_ptr<BackendType> backend_;
  mutable std::mutex mutex_;
}; // class CommonDenseMatrix

/**
 * \brief A sparse matrix implementation of the MatrixInterface with row major memory layout.
 */
template <class ScalarImp = double>
class CommonSparseMatrix : public MatrixInterface<internal::CommonSparseMatrixTraits<ScalarImp>, ScalarImp>
{
  typedef CommonSparseMatrix<ScalarImp> ThisType;
  typedef MatrixInterface<internal::CommonSparseMatrixTraits<ScalarImp>, ScalarImp> MatrixInterfaceType;

public:
  typedef internal::CommonSparseMatrixTraits<ScalarImp> Traits;
  typedef typename Traits::EntriesVectorType EntriesVectorType;
  typedef typename Traits::ScalarType ScalarType;
  typedef typename Traits::RealType RealType;
  typedef std::vector<size_t> IndexVectorType;

  /**
  * \brief This is the constructor of interest which creates a sparse matrix.
  */
  CommonSparseMatrix(const size_t rr, const size_t cc, const SparsityPatternDefault& patt)
    : num_rows_(rr)
    , num_cols_(cc)
    , entries_(std::make_shared<EntriesVectorType>())
    , row_pointers_(std::make_shared<IndexVectorType>(num_rows_ + 1))
    , column_indices_(std::make_shared<IndexVectorType>())
  {
    if (num_rows_ > 0 && num_cols_ > 0) {
      if (size_t(patt.size()) != num_rows_)
        DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                   "The size of the pattern (" << patt.size() << ") does not match the number of rows of this ("
                                               << num_rows_
                                               << ")!");
      for (size_t row = 0; row < num_rows_; ++row) {
        const auto& columns = patt.inner(row);
        const auto num_nonzero_entries_in_row = columns.size();
        row_pointers_->operator[](row + 1) = row_pointers_->operator[](row) + num_nonzero_entries_in_row;
        for (size_t col = 0; col < num_nonzero_entries_in_row; ++col) {
#ifndef NDEBUG
          if (col >= num_cols_)
            DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                       "The size of row " << row << " of the pattern does not match the number of columns of this ("
                                          << num_cols_
                                          << ")!");
#endif // NDEBUG
          column_indices_->push_back(columns[col]);
        }
        entries_->resize(column_indices_->size());
      }
    }
  }

  CommonSparseMatrix(const size_t rr = 0, const size_t cc = 0, const ScalarType& value = ScalarType(0))
    : num_rows_(rr)
    , num_cols_(cc)
    , entries_(std::make_shared<EntriesVectorType>(num_rows_ * num_cols_, value))
    , row_pointers_(std::make_shared<IndexVectorType>(num_rows_ + 1))
    , column_indices_(std::make_shared<IndexVectorType>())
  {
    IndexVectorType row_column_indices(num_cols_);
    for (size_t col = 0; col < num_cols_; ++col) {
      row_column_indices[col] = col;
    }
    for (size_t row = 0; row < num_rows_; ++row) {
      row_pointers_->operator[](row + 1) = (row + 1) * num_cols_;
      column_indices_->insert(column_indices_->end(), row_column_indices.begin(), row_column_indices.end());
    }
  }

  explicit CommonSparseMatrix(const int rr, const int cc = 0)
    : CommonSparseMatrix(size_t(rr), size_t(cc), ScalarType(0))
  {
    assert(rr >= 0 && cc >= 0);
  }

  CommonSparseMatrix(const ThisType& other)
    : num_rows_(other.num_rows_)
    , num_cols_(other.num_cols_)
    , entries_(other.entries_)
    , row_pointers_(other.row_pointers_)
    , column_indices_(other.column_indices_)
  {
  }

  template <class OtherMatrixType>
  explicit CommonSparseMatrix(
      const OtherMatrixType& mat,
      const typename std::enable_if<Common::MatrixAbstraction<OtherMatrixType>::is_matrix, bool>::type prune = false,
      const typename Common::FloatCmp::DefaultEpsilon<ScalarType>::Type eps =
          Common::FloatCmp::DefaultEpsilon<ScalarType>::value())
    : num_rows_(Common::MatrixAbstraction<OtherMatrixType>::rows(mat))
    , num_cols_(Common::MatrixAbstraction<OtherMatrixType>::cols(mat))
    , entries_(std::make_shared<EntriesVectorType>())
    , row_pointers_(std::make_shared<IndexVectorType>(num_rows_ + 1))
    , column_indices_(std::make_shared<IndexVectorType>())
  {
    for (size_t rr = 0; rr < num_rows_; ++rr) {
      size_t num_nonzero_entries_in_row = 0;
      for (size_t cc = 0; cc < num_cols_; ++cc) {
        const auto& value = Common::MatrixAbstraction<OtherMatrixType>::get_entry(mat, rr, cc);
        if (!prune || XT::Common::FloatCmp::ne(value, ScalarType(0), eps)) {
          ++num_nonzero_entries_in_row;
          entries_->push_back(value);
          column_indices_->push_back(cc);
        }
      }
      row_pointers_->operator[](rr + 1) = row_pointers_->operator[](rr) + num_nonzero_entries_in_row;
    }
  } // CommonSparseMatrix(...)

  template <int ROWS, int COLS>
  explicit operator Dune::FieldMatrix<ScalarType, ROWS, COLS>() const
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    assert(ROWS == num_rows_ && COLS == num_cols_);
    Dune::FieldMatrix<ScalarType, ROWS, COLS> ret(ScalarType(0));
    for (size_t rr = 0; rr < ROWS; ++rr)
      for (size_t kk = row_pointers_->operator[](rr); kk < row_pointers_->operator[](rr + 1); ++kk)
        ret[rr][column_indices_->operator[](kk)] = entries_->operator[](kk);
    return ret;
  }

  ThisType& operator=(const ThisType& other)
  {
    entries_ = other.entries_;
    row_pointers_ = other.row_pointers_;
    column_indices_ = other.column_indices_;
    num_rows_ = other.num_rows_;
    num_cols_ = other.num_cols_;
    return *this;
  }

  /// \}
  /// \name Required by ContainerInterface.
  /// \{

  inline ThisType copy() const
  {
    ThisType ret(*this);
    return ret;
  }

  inline void scal(const ScalarType& alpha)
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    ensure_uniqueness();
    std::transform(
        entries_->begin(), entries_->end(), entries_->begin(), std::bind1st(std::multiplies<ScalarType>(), alpha));
  }

  inline void axpy(const ScalarType& alpha, const ThisType& xx)
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    ensure_uniqueness();
    assert(has_equal_shape(xx));
    const auto& xx_entries = *xx.entries_;
    for (size_t ii = 0; ii < entries_->size(); ++ii)
      entries_->operator[](ii) += alpha * xx_entries[ii];
  }

  inline bool has_equal_shape(const ThisType& other) const
  {
    return (rows() == other.rows()) && (cols() == other.cols());
  }

  /// \}
  /// \name Required by MatrixInterface.
  /// \{

  inline size_t rows() const
  {
    return num_rows_;
  }

  inline size_t cols() const
  {
    return num_cols_;
  }

  template <class XX, class YY>
  inline void mv(const XX& xx, YY& yy) const
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    std::fill(yy.begin(), yy.end(), ScalarType(0));
    for (size_t rr = 0; rr < num_rows_; ++rr)
      for (size_t kk = row_pointers_->operator[](rr); kk < row_pointers_->operator[](rr + 1); ++kk)
        yy[rr] += entries_->operator[](kk) * xx[column_indices_->operator[](kk)];
  }

  inline void add_to_entry(const size_t rr, const size_t cc, const ScalarType& value)
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    ensure_uniqueness();
    entries_->operator[](get_entry_index(rr, cc)) += value;
  }

  inline ScalarType get_entry(const size_t rr, const size_t cc) const
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    const size_t index = get_entry_index(rr, cc, false);
    return index == size_t(-1) ? ScalarType(0) : entries_->operator[](index);
  }

  inline ScalarType& get_entry_ref(const size_t rr, const size_t cc)
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    const size_t index = get_entry_index(rr, cc, false);
    if (index == size_t(-1))
      DUNE_THROW(Dune::RangeError, "Entry not in matrix pattern!");
    return entries_->operator[](index);
  }

  inline void set_entry(const size_t rr, const size_t cc, const ScalarType value)
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    ensure_uniqueness();
    entries_->operator[](get_entry_index(rr, cc)) = value;
  }

  inline void clear_row(const size_t rr)
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    ensure_uniqueness();
    std::fill(entries_->begin() + row_pointers_->operator[](rr),
              entries_->begin() + row_pointers_->operator[](rr + 1),
              ScalarType(0));
  }

  inline void clear_col(const size_t cc)
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    ensure_uniqueness();
    for (size_t kk = 0; kk < entries_->size(); ++kk) {
      if (column_indices_->operator[](kk) == cc)
        entries_->operator[](kk) = ScalarType(0);
    }
  }

  inline void unit_row(const size_t rr)
  {
    ensure_uniqueness();
    clear_row(rr);
    set_entry(rr, rr, ScalarType(1));
  }

  inline void unit_col(const size_t cc)
  {
    ensure_uniqueness();
    clear_col(cc);
    set_entry(cc, cc, ScalarType(1));
  }

  bool valid() const
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    // iterate over non-zero entries
    for (const auto& entry : *entries_)
      if (XT::Common::isnan(std::real(entry)) || XT::Common::isnan(std::imag(entry))
          || XT::Common::isinf(std::abs(entry)))
        return false;
    return true;
  }

  virtual size_t non_zeros() const override final
  {
    return entries_->size();
  }

  virtual SparsityPatternDefault pattern(const bool prune = false,
                                         const typename Common::FloatCmp::DefaultEpsilon<ScalarType>::Type eps =
                                             Common::FloatCmp::DefaultEpsilon<ScalarType>::value()) const override
  {
    SparsityPatternDefault ret(num_rows_);
    for (size_t rr = 0; rr < num_rows_; ++rr) {
      for (size_t kk = row_pointers_->operator[](rr + 1); kk < row_pointers_->operator[](rr + 1); ++kk) {
        const auto& val = entries_->operator[](kk);
        if (!prune || Common::FloatCmp::ne(val, ScalarType(0), eps))
          ret.insert(rr, column_indices_->operator[](kk));
      }
    }
    return ret;
  } // ... pattern(...)

  /// \}

private:
  size_t get_entry_index(const size_t rr, const size_t cc, const bool throw_if_not_in_pattern = true) const
  {
    const auto& row_offset = row_pointers_->operator[](rr);
    const auto& next_row_offset = row_pointers_->operator[](rr + 1);
    const auto column_indices_it = column_indices_->begin() + row_offset;
    const auto column_indices_it_end = column_indices_->begin() + next_row_offset;
    const auto entry_it = std::find(column_indices_it, column_indices_it_end, cc);
    if (entry_it != column_indices_it_end)
      return row_offset + std::distance(column_indices_it, entry_it);
    if (throw_if_not_in_pattern)
      DUNE_THROW(Common::Exceptions::index_out_of_range, "Entry is not in the sparsity pattern!");
    return size_t(-1);
  }

protected:
  inline void ensure_uniqueness()
  {
    if (!entries_.unique())
      entries_ = std::make_shared<EntriesVectorType>(*entries_);
  } // ... ensure_uniqueness(...)

private:
  size_t num_rows_, num_cols_;
  std::shared_ptr<EntriesVectorType> entries_;
  std::shared_ptr<IndexVectorType> row_pointers_;
  std::shared_ptr<IndexVectorType> column_indices_;
  mutable std::mutex mutex_;
}; // class CommonSparseMatrix

} // namespace LA
namespace Common {

template <class T>
struct VectorAbstraction<LA::CommonDenseVector<T>>
    : public LA::internal::VectorAbstractionBase<LA::CommonDenseVector<T>>
{
};

template <class T>
struct MatrixAbstraction<LA::CommonDenseMatrix<T>>
    : public LA::internal::MatrixAbstractionBase<LA::CommonDenseMatrix<T>>
{
};

template <class T>
struct MatrixAbstraction<LA::CommonSparseMatrix<T>>
    : public LA::internal::MatrixAbstractionBase<LA::CommonSparseMatrix<T>>
{
};

} // namespace Common
} // namespace XT
} // namespace Dune


// begin: this is what we need for the lib
#if DUNE_XT_WITH_PYTHON_BINDINGS


extern template class Dune::XT::LA::CommonDenseVector<double>;
extern template class Dune::XT::LA::CommonDenseMatrix<double>;
extern template class Dune::XT::LA::CommonSparseMatrix<double>;
// extern template void Dune::XT::LA::CommonSparseMatrix<double>::mv(const Dune::XT::LA::CommonDenseVector<double>&,
//                                                                  Dune::XT::LA::CommonDenseVector<double>) const;


#endif // DUNE_XT_WITH_PYTHON_BINDINGS
// end: this is what we need for the lib


#endif // DUNE_XT_LA_CONTAINER_COMMON_HH
