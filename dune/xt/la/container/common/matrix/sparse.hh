// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Barbara Verfürth (2015)
//   Felix Schindler  (2013 - 2017)
//   Rene Milk        (2014 - 2016, 2018)
//   TiKeil           (2018)
//   Tobias Leibner   (2014, 2016 - 2018)

#ifndef DUNE_XT_LA_CONTAINER_COMMON_MATRIX_SPARSE_HH
#define DUNE_XT_LA_CONTAINER_COMMON_MATRIX_SPARSE_HH

#include <dune/xt/common/matrix.hh>

#include <dune/xt/la/container/interfaces.hh>
#include <dune/xt/la/container/pattern.hh>

#include "../vector/sparse.hh"
#include "dense.hh"

namespace Dune {
namespace XT {
namespace LA {


// forwards
template <class ScalarImp, Common::StorageLayout layout>
class CommonSparseMatrix;

// forwards
template <class DenseMatrixImp, class SparseMatrixImp>
class CommonSparseOrDenseMatrix;


namespace internal {


template <class ScalarImp, Common::StorageLayout layout>
struct CommonSparseMatrixTraits : public MatrixTraitsBase<ScalarImp,
                                                          CommonSparseMatrix<ScalarImp, layout>,
                                                          void,
                                                          Backends::common_sparse,
                                                          Backends::common_dense,
                                                          true>
{
  using EntriesVectorType = std::vector<ScalarImp>;
  using IndexVectorType = std::vector<size_t>;
  using EpsType = typename Common::FloatCmp::DefaultEpsilon<ScalarImp>::Type;
};


template <class DenseMatrixImp, class SparseMatrixImp>
struct CommonSparseOrDenseMatrixTraits
    : public MatrixTraitsBase<typename DenseMatrixImp::ScalarType,
                              CommonSparseOrDenseMatrix<DenseMatrixImp, SparseMatrixImp>,
                              void,
                              Backends::common_dense,
                              Backends::common_dense,
                              false>
{
  using DenseMatrixType = DenseMatrixImp;
  using SparseMatrixType = SparseMatrixImp;
  using EpsType = typename Common::FloatCmp::DefaultEpsilon<typename DenseMatrixImp::ScalarType>::Type;
};


} // namespace internal


/**
 * \brief A sparse matrix implementation of the MatrixInterface with row major memory layout.
 */
template <class ScalarImp = double, Common::StorageLayout layout = Common::StorageLayout::csr>
class CommonSparseMatrix : public MatrixInterface<internal::CommonSparseMatrixTraits<ScalarImp, layout>, ScalarImp>
{
  using ThisType = CommonSparseMatrix;
  using InterfaceType = MatrixInterface<internal::CommonSparseMatrixTraits<ScalarImp, layout>, ScalarImp>;

public:
  using typename InterfaceType::RealType;
  using typename InterfaceType::ScalarType;
  using typename InterfaceType::Traits;
  using EntriesVectorType = typename Traits::EntriesVectorType;
  using EpsType = typename Traits::EpsType;
  using IndexVectorType = typename Traits::IndexVectorType;

private:
  using MutexesType = typename Traits::MutexesType;

public:
  /**
  * \brief This is the constructor of interest which creates a sparse matrix.
  */
  CommonSparseMatrix(const size_t rr,
                     const size_t cc,
                     const SparsityPatternDefault& patt,
                     const size_t num_mutexes = 1,
                     const EpsType eps = Common::FloatCmp::DefaultEpsilon<ScalarType>::value() / 1000.)
    : num_rows_(rr)
    , num_cols_(cc)
    , entries_(std::make_shared<EntriesVectorType>())
    , row_pointers_(std::make_shared<IndexVectorType>(num_rows_ + 1, 0))
    , column_indices_(std::make_shared<IndexVectorType>())
    , mutexes_(std::make_shared<std::vector<std::mutex>>(num_mutexes))
    , eps_(eps)
    , unshareable_(false)
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
        assert(columns.size() <= num_cols_);
        row_pointers_->operator[](row + 1) = row_pointers_->operator[](row) + num_nonzero_entries_in_row;
        for (size_t kk = 0; kk < num_nonzero_entries_in_row; ++kk) {
#ifndef NDEBUG
          if (columns[kk] >= num_cols_)
            DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                       "The size of row " << row << " of the pattern does not match the number of columns of this ("
                                          << num_cols_
                                          << ")!");
#endif // NDEBUG
          column_indices_->push_back(columns[kk]);
        } // kk
        entries_->resize(column_indices_->size());
      } // row
    }
  } // CommonSparseMatrix(rr, cc, patt, num_mutexes)

  CommonSparseMatrix(const size_t rr = 0,
                     const size_t cc = 0,
                     const ScalarType& value = ScalarType(0),
                     const size_t num_mutexes = 1,
                     const EpsType eps = Common::FloatCmp::DefaultEpsilon<ScalarType>::value() / 1000.)
    : num_rows_(rr)
    , num_cols_(cc)
    , entries_(std::make_shared<EntriesVectorType>(is_zero(value, eps) ? 0 : num_rows_ * num_cols_, value))
    , row_pointers_(std::make_shared<IndexVectorType>(num_rows_ + 1, 0))
    , column_indices_(std::make_shared<IndexVectorType>())
    , mutexes_(std::make_shared<std::vector<std::mutex>>(num_mutexes))
    , eps_(eps)
    , unshareable_(false)
  {
    if (!is_zero(value, eps)) {
      IndexVectorType row_column_indices(num_cols_);
      for (size_t col = 0; col < num_cols_; ++col)
        row_column_indices[col] = col;
      for (size_t row = 0; row < num_rows_; ++row) {
        row_pointers_->operator[](row + 1) = (row + 1) * num_cols_;
        column_indices_->insert(column_indices_->end(), row_column_indices.begin(), row_column_indices.end());
      }
    }
  }

  CommonSparseMatrix(const ThisType& other)
    : num_rows_(other.num_rows_)
    , num_cols_(other.num_cols_)
    , entries_(other.unshareable_ ? std::make_shared<EntriesVectorType>(*other.entries_) : other.entries_)
    , row_pointers_(other.row_pointers_)
    , column_indices_(other.column_indices_)
    , mutexes_(other.unshareable_ ? std::make_shared<std::vector<std::mutex>>(other.mutexes_->size()) : other.mutexes_)
    , eps_(other.eps_)
    , unshareable_(false)
  {
  }

  template <class OtherMatrixType>
  explicit CommonSparseMatrix(
      const OtherMatrixType& mat,
      const typename std::enable_if<Common::MatrixAbstraction<OtherMatrixType>::is_matrix, bool>::type prune = false,
      const EpsType eps = Common::FloatCmp::DefaultEpsilon<ScalarType>::value() / 1000.,
      const size_t num_mutexes = 1)
    : num_rows_(Common::MatrixAbstraction<OtherMatrixType>::rows(mat))
    , num_cols_(Common::MatrixAbstraction<OtherMatrixType>::cols(mat))
    , entries_(std::make_shared<EntriesVectorType>())
    , row_pointers_(std::make_shared<IndexVectorType>(num_rows_ + 1))
    , column_indices_(std::make_shared<IndexVectorType>())
    , mutexes_(std::make_shared<std::vector<std::mutex>>(num_mutexes))
    , eps_(eps)
    , unshareable_(false)
  {
    size_t index = 0;
    for (size_t rr = 0; rr < num_rows_; ++rr) {
      (*row_pointers_)[rr] = index;
      for (size_t cc = 0; cc < num_cols_; ++cc) {
        const auto& value = Common::MatrixAbstraction<OtherMatrixType>::get_entry(mat, rr, cc);
        if (!prune || XT::Common::FloatCmp::ne(
                          value, ScalarType(0), 0., eps / Common::MatrixAbstraction<OtherMatrixType>::cols(mat))) {
          entries_->push_back(value);
          column_indices_->push_back(cc);
          ++index;
        }
      } // cc
    } // rr
    (*row_pointers_)[num_rows_] = column_indices_->size();
  } // CommonSparseMatrix(...)

  template <class OtherMatrixType>
  explicit CommonSparseMatrix(
      const OtherMatrixType& mat,
      const typename std::enable_if<Common::MatrixAbstraction<OtherMatrixType>::is_matrix, bool>::type prune,
      const size_t num_mutexes)
    : CommonSparseMatrix(mat,
                         prune,
                         Common::FloatCmp::DefaultEpsilon<ScalarType>::value()
                             / (1000. * Common::MatrixAbstraction<OtherMatrixType>::cols(mat)),
                         num_mutexes)
  {
  } // CommonSparseMatrix(...)

  template <class DenseMatrixImp>
  void copy_to_densematrix(DenseMatrixImp& ret) const
  {
    for (size_t rr = 0; rr < num_rows_; ++rr)
      for (size_t kk = row_pointers_->operator[](rr); kk < row_pointers_->operator[](rr + 1); ++kk)
        ret[rr][column_indices_->operator[](kk)] = entries_->operator[](kk);
  }

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      num_rows_ = other.num_rows_;
      num_cols_ = other.num_cols_;
      entries_ = other.unshareable_ ? std::make_shared<EntriesVectorType>(*other.entries_) : other.entries_;
      row_pointers_ = other.row_pointers_;
      column_indices_ = other.column_indices_;
      mutexes_ =
          other.unshareable_ ? std::make_shared<std::vector<std::mutex>>(other.mutexes_->size()) : other.mutexes_;
      unshareable_ = false;
    }
    return *this;
  }

  template <class OtherMatrixImp>
  typename std::enable_if_t<XT::Common::MatrixAbstraction<OtherMatrixImp>::is_matrix, ThisType>&
  assign(const OtherMatrixImp& other, const SparsityPatternDefault& pattern)
  {
    clear();
    using MatAbstrType = XT::Common::MatrixAbstraction<OtherMatrixImp>;
    num_rows_ = MatAbstrType::rows(other);
    num_cols_ = MatAbstrType::cols(other);
    for (size_t rr = 0; rr < num_rows_; ++rr) {
      for (const auto& cc : pattern.inner(rr)) {
        if (XT::Common::FloatCmp::ne(MatAbstrType::get_entry(other, rr, cc), 0., 0., eps_ / num_cols_)) {
          entries_->push_back(MatAbstrType::get_entry(other, rr, cc));
          column_indices_->push_back(cc);
        }
      } // cc
      (*row_pointers_)[rr + 1] = column_indices_->size();
    } // rr
    return *this;
  }

  template <class OtherMatrixImp>
  typename std::enable_if_t<XT::Common::MatrixAbstraction<OtherMatrixImp>::is_matrix, ThisType>&
  operator=(const OtherMatrixImp& other)
  {
    return assign(other, dense_pattern(num_rows_, num_cols_));
  }

  void deep_copy(const ThisType& other)
  {
    ensure_uniqueness();
    num_rows_ = other.num_rows_;
    num_cols_ = other.num_cols_;
    *entries_ = *other.entries_;
    *row_pointers_ = *other.row_pointers_;
    *column_indices_ = *other.column_indices_;
  }

  void clear()
  {
    ensure_uniqueness();
    entries_->clear();
    std::fill(row_pointers_->begin(), row_pointers_->end(), 0);
    column_indices_->clear();
  }

  /// \name Required by ContainerInterface.
  /// \{
  inline ThisType copy() const
  {
    ThisType ret(*this);
    ret.entries_ = std::make_shared<EntriesVectorType>(*entries_);
    ret.row_pointers_ = std::make_shared<IndexVectorType>(*row_pointers_);
    ret.column_indices_ = std::make_shared<IndexVectorType>(*column_indices_);
    ret.mutexes_ = std::make_shared<std::vector<std::mutex>>(mutexes_->size());
    return ret;
  }

  inline void scal(const ScalarType& alpha)
  {
    ensure_uniqueness();
    const internal::VectorLockGuard DUNE_UNUSED(guard)(*mutexes_);
    std::transform(
        entries_->begin(), entries_->end(), entries_->begin(), std::bind1st(std::multiplies<ScalarType>(), alpha));
  }

  inline void axpy(const ScalarType& alpha, const ThisType& xx)
  {
    ensure_uniqueness();
    const internal::VectorLockGuard DUNE_UNUSED(guard)(*mutexes_);
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

  //! Matrix-Vector multiplication for arbitrary vectors that support operator[]
  template <class XX, class YY>
  inline std::enable_if_t<XT::Common::VectorAbstraction<XX>::is_vector && XT::Common::VectorAbstraction<YY>::is_vector,
                          void>
  mv(const XX& xx, YY& yy) const
  {
    std::fill(yy.begin(), yy.end(), ScalarType(0));
    const auto& entries = *entries_;
    const auto& row_pointers = *row_pointers_;
    const auto& column_indices = *column_indices_;
    for (size_t rr = 0; rr < num_rows_; ++rr) {
      auto& yy_rr = yy[rr];
      const size_t end = row_pointers[rr + 1];
      for (size_t kk = row_pointers[rr]; kk < end; ++kk)
        yy_rr += entries[kk] * xx[column_indices[kk]];
    }
  }

  //! TransposedMatrix-Vector multiplication for arbitrary vectors that support operator[]
  template <class XX, class YY>
  inline std::enable_if_t<XT::Common::VectorAbstraction<XX>::is_vector && XT::Common::VectorAbstraction<YY>::is_vector,
                          void>
  mtv(const XX& xx, YY& yy) const
  {
    std::fill(yy.begin(), yy.end(), ScalarType(0));
    const auto& entries = *entries_;
    const auto& row_pointers = *row_pointers_;
    const auto& column_indices = *column_indices_;
    for (size_t rr = 0; rr < num_rows_; ++rr) {
      const size_t end = row_pointers[rr + 1];
      for (size_t kk = row_pointers[rr]; kk < end; ++kk)
        yy[column_indices[kk]] += entries[kk] * xx[rr];
    }
  }

  inline void add_to_entry(const size_t rr, const size_t cc, const ScalarType& value)
  {
    ensure_uniqueness();
    internal::LockGuard DUNE_UNUSED(lock)(*mutexes_, rr, rows());
    entries_->operator[](get_entry_index(rr, cc)) += value;
  }

  inline ScalarType get_entry(const size_t rr, const size_t cc) const
  {
    const size_t index = get_entry_index(rr, cc, false);
    return index == size_t(-1) ? ScalarType(0) : entries_->operator[](index);
  }

  inline void set_entry(const size_t rr, const size_t cc, const ScalarType value)
  {
    ensure_uniqueness();
    entries_->operator[](get_entry_index(rr, cc)) = value;
  }

  inline void start_row()
  {
    if (row_pointers_->empty())
      row_pointers_->push_back(0);
  }

  inline void end_row()
  {
    row_pointers_->push_back(column_indices_->size());
  }

  inline void push_entry(const size_t cc, const ScalarType value)
  {
    ensure_uniqueness();
    entries_->push_back(value);
    column_indices_->push_back(cc);
  }

  inline void clear_row(const size_t rr)
  {
    ensure_uniqueness();
    std::fill(entries_->begin() + row_pointers_->operator[](rr),
              entries_->begin() + row_pointers_->operator[](rr + 1),
              ScalarType(0));
  }

  inline void clear_col(const size_t cc)
  {
    ensure_uniqueness();
    for (size_t kk = 0; kk < entries_->size(); ++kk) {
      if (column_indices_->operator[](kk) == cc)
        entries_->operator[](kk) = ScalarType(0);
    }
  }

  inline void unit_row(const size_t rr)
  {
    clear_row(rr);
    set_entry(rr, rr, ScalarType(1));
  }

  inline void unit_col(const size_t cc)
  {
    clear_col(cc);
    set_entry(cc, cc, ScalarType(1));
  }

  bool valid() const
  {
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
                                         const EpsType eps = Common::FloatCmp::DefaultEpsilon<ScalarType>::value()
                                                             / 1000.) const override
  {
    SparsityPatternDefault ret(num_rows_);
    for (size_t rr = 0; rr < num_rows_; ++rr) {
      for (size_t kk = row_pointers_->operator[](rr); kk < row_pointers_->operator[](rr + 1); ++kk) {
        if (!prune || !is_zero(entries_->operator[](kk), eps))
          ret.insert(rr, column_indices_->operator[](kk));
      }
    }
    return ret;
  } // ... pattern(...)

  /// \}

  template <class MatrixType>
  void rightmultiply(const MatrixType& other)
  {
    ensure_uniqueness();
    auto new_entries = std::make_shared<EntriesVectorType>();
    auto new_row_pointers = std::make_shared<IndexVectorType>(num_rows_ + 1);
    auto new_column_indices = std::make_shared<IndexVectorType>();
    new_entries->reserve(entries_->size());
    new_column_indices->reserve(column_indices_->size());
    ScalarType new_entry(0);
    size_t index = 0;
    for (size_t rr = 0; rr < num_rows_; ++rr) {
      (*new_row_pointers)[rr] = index;
      for (size_t cc = 0; cc < num_cols_; ++cc) {
        new_entry = 0;
        for (size_t kk = (*row_pointers_)[rr]; kk < (*row_pointers_)[rr + 1]; ++kk)
          new_entry += (*entries_)[kk] * other[(*column_indices_)[kk]][cc];
        if (XT::Common::FloatCmp::ne(new_entry, 0., 0., eps_ / num_cols_)) {
          new_entries->push_back(new_entry);
          new_column_indices->push_back(cc);
          ++index;
        }
      } // cc
    } // rr
    (*new_row_pointers)[num_rows_] = new_column_indices->size();
    entries_ = new_entries;
    row_pointers_ = new_row_pointers;
    column_indices_ = new_column_indices;
  } // void rightmultiply(...)

  using InterfaceType::operator+;
  using InterfaceType::operator-;
  using InterfaceType::operator+=;
  using InterfaceType::operator-=;

  ScalarType* entries()
  {
    ensure_uniqueness();
    unshareable_ = true;
    return entries_->data();
  }

  const ScalarType* entries() const
  {
    ensure_uniqueness();
    unshareable_ = true;
    return entries_->data();
  }

  size_t* outer_index_ptr()
  {
    ensure_uniqueness();
    unshareable_ = true;
    return row_pointers_->data();
  }

  const size_t* outer_index_ptr() const
  {
    ensure_uniqueness();
    unshareable_ = true;
    return row_pointers_->data();
  }

  size_t* inner_index_ptr()
  {
    ensure_uniqueness();
    unshareable_ = true;
    return column_indices_->data();
  }

  const size_t* inner_index_ptr() const
  {
    ensure_uniqueness();
    unshareable_ = true;
    return column_indices_->data();
  }

  const EpsType& eps()
  {
    return eps_;
  }

  inline void ensure_uniqueness() const
  {
    if (!entries_.unique()) {
      assert(!unshareable_);
      const internal::VectorLockGuard DUNE_UNUSED(guard)(*mutexes_);
      if (!entries_.unique()) {
        entries_ = std::make_shared<EntriesVectorType>(*entries_);
        row_pointers_ = std::make_shared<IndexVectorType>(*row_pointers_);
        column_indices_ = std::make_shared<IndexVectorType>(*column_indices_);
        mutexes_ = std::make_shared<std::vector<std::mutex>>(mutexes_->size());
      }
    }
  } // ... ensure_uniqueness(...)

private:
  size_t get_entry_index(const size_t rr, const size_t cc, const bool throw_if_not_in_pattern = true) const
  {
    const auto& row_offset = row_pointers_->operator[](rr);
    const auto& next_row_offset = row_pointers_->operator[](rr + 1);
    const auto column_indices_it = column_indices_->begin() + row_offset;
    const auto column_indices_it_end = column_indices_->begin() + next_row_offset;
    const auto entry_it = std::lower_bound(column_indices_it, column_indices_it_end, cc);
    if (entry_it != column_indices_it_end && *entry_it == cc)
      return row_offset + std::distance(column_indices_it, entry_it);
    if (throw_if_not_in_pattern)
      DUNE_THROW(Common::Exceptions::index_out_of_range, "Entry is not in the sparsity pattern!");
    return size_t(-1);
  }

  bool is_zero(const ScalarType val, const ScalarType eps) const
  {
    const auto factor = (num_cols_ == 0 ? 1. : ScalarType(num_cols_));
    const ScalarType tol = eps / factor;
    return XT::Common::FloatCmp::eq(val, ScalarType(0.), 0., tol);
  }

  size_t num_rows_, num_cols_;
  mutable std::shared_ptr<EntriesVectorType> entries_;
  mutable std::shared_ptr<IndexVectorType> row_pointers_;
  mutable std::shared_ptr<IndexVectorType> column_indices_;
  mutable std::shared_ptr<std::vector<std::mutex>> mutexes_;
  mutable EpsType eps_;
  mutable bool unshareable_;
}; // class CommonSparseMatrix

/**
 * \brief A sparse matrix implementation of the MatrixInterface with column major memory layout.
 */
template <class ScalarImp>
class CommonSparseMatrix<ScalarImp, Common::StorageLayout::csc>
    : public MatrixInterface<internal::CommonSparseMatrixTraits<ScalarImp, Common::StorageLayout::csc>, ScalarImp>
{
  using ThisType = CommonSparseMatrix;
  using InterfaceType =
      MatrixInterface<internal::CommonSparseMatrixTraits<ScalarImp, Common::StorageLayout::csc>, ScalarImp>;

public:
  using typename InterfaceType::RealType;
  using typename InterfaceType::ScalarType;
  using typename InterfaceType::Traits;
  using EntriesVectorType = typename Traits::EntriesVectorType;
  using EpsType = typename Traits::EpsType;
  using IndexVectorType = typename Traits::IndexVectorType;

private:
  using MutexesType = typename Traits::MutexesType;

public:
  /**
  * \brief This is the constructor of interest which creates a sparse matrix.
  */
  CommonSparseMatrix(const size_t rr,
                     const size_t cc,
                     const SparsityPatternDefault& patt,
                     const size_t num_mutexes = 1,
                     const EpsType eps = Common::FloatCmp::DefaultEpsilon<ScalarType>::value() / 1000.)
    : num_rows_(rr)
    , num_cols_(cc)
    , entries_(std::make_shared<EntriesVectorType>())
    , column_pointers_(std::make_shared<IndexVectorType>(num_cols_ + 1, 0))
    , row_indices_(std::make_shared<IndexVectorType>())
    , mutexes_(std::make_shared<std::vector<std::mutex>>(num_mutexes))
    , eps_(eps)
    , unshareable_(false)
  {
    if (num_rows_ > 0 && num_cols_ > 0) {
      if (size_t(patt.size()) != num_rows_)
        DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                   "The size of the pattern (" << patt.size() << ") does not match the number of rows of this ("
                                               << num_rows_
                                               << ")!");
      for (size_t col = 0; col < num_cols_; ++col) {
        for (size_t row = 0; row < num_rows_; ++row) {
          const auto& column_indices = patt.inner(row);
          if (std::find(column_indices.begin(), column_indices.end(), col) != column_indices.end())
            row_indices_->push_back(row);
        } // row
        (*column_pointers_)[col + 1] = row_indices_->size();
      } // col
      entries_->resize(row_indices_->size());
    }
  } // CommonSparseMatrix(rr, cc, patt, num_mutexes)

  CommonSparseMatrix(const size_t rr = 0,
                     const size_t cc = 0,
                     const ScalarType& value = ScalarType(0),
                     const size_t num_mutexes = 1,
                     const EpsType eps = Common::FloatCmp::DefaultEpsilon<ScalarType>::value() / 1000.)
    : num_rows_(rr)
    , num_cols_(cc)
    , entries_(std::make_shared<EntriesVectorType>(is_zero(value, eps) ? 0 : num_rows_ * num_cols_, value))
    , column_pointers_(std::make_shared<IndexVectorType>(num_cols_ + 1, 0))
    , row_indices_(std::make_shared<IndexVectorType>())
    , mutexes_(std::make_shared<std::vector<std::mutex>>(num_mutexes))
    , eps_(eps)
    , unshareable_(false)
  {
    if (!is_zero(value, eps)) {
      IndexVectorType column_row_indices(num_rows_);
      for (size_t row = 0; row < num_rows_; ++row)
        column_row_indices[row] = row;
      for (size_t col = 0; col < num_cols_; ++col) {
        column_pointers_->operator[](col + 1) = (col + 1) * num_rows_;
        row_indices_->insert(row_indices_->end(), column_row_indices.begin(), column_row_indices.end());
      }
    }
  }

  CommonSparseMatrix(const size_t rr, const size_t cc, const size_t num_mutexes)
    : num_rows_(rr)
    , num_cols_(cc)
    , entries_(std::make_shared<EntriesVectorType>())
    , column_pointers_(std::make_shared<IndexVectorType>(num_cols_ + 1))
    , row_indices_(std::make_shared<IndexVectorType>())
    , mutexes_(std::make_shared<std::vector<std::mutex>>(num_mutexes))
    , unshareable_(false)
  {
  }

  CommonSparseMatrix(const ThisType& other)
    : num_rows_(other.num_rows_)
    , num_cols_(other.num_cols_)
    , entries_(other.unshareable_ ? std::make_shared<EntriesVectorType>(*other.entries_) : other.entries_)
    , column_pointers_(other.unshareable_ ? std::make_shared<IndexVectorType>(*other.column_pointers_)
                                          : other.column_pointers_)
    , row_indices_(other.unshareable_ ? std::make_shared<IndexVectorType>(*other.row_indices_) : other.row_indices_)
    , mutexes_(other.unshareable_ ? std::make_shared<std::vector<std::mutex>>(other.mutexes_->size()) : other.mutexes_)
    , unshareable_(false)
  {
  }

  template <class OtherMatrixType>
  explicit CommonSparseMatrix(
      const OtherMatrixType& mat,
      const typename std::enable_if<Common::MatrixAbstraction<OtherMatrixType>::is_matrix, bool>::type prune = false,
      const EpsType eps = Common::FloatCmp::DefaultEpsilon<ScalarType>::value() / 1000.,
      const size_t num_mutexes = 1)
    : num_rows_(Common::MatrixAbstraction<OtherMatrixType>::rows(mat))
    , num_cols_(Common::MatrixAbstraction<OtherMatrixType>::cols(mat))
    , entries_(std::make_shared<EntriesVectorType>())
    , column_pointers_(std::make_shared<IndexVectorType>(num_cols_ + 1, 0))
    , row_indices_(std::make_shared<IndexVectorType>())
    , mutexes_(std::make_shared<std::vector<std::mutex>>(num_mutexes))
    , eps_(eps)
    , unshareable_(false)
  {
    entries_->reserve(0.1 * num_rows_ * num_cols_);
    row_indices_->reserve(0.1 * num_rows_ * num_cols_);
    size_t index = 0;
    for (size_t cc = 0; cc < num_cols_; ++cc) {
      for (size_t rr = 0; rr < num_rows_; ++rr) {
        const auto& value = Common::MatrixAbstraction<OtherMatrixType>::get_entry(mat, rr, cc);
        if (!prune || XT::Common::FloatCmp::ne(
                          value, ScalarType(0), 0., eps / Common::MatrixAbstraction<OtherMatrixType>::cols(mat))) {
          entries_->push_back(value);
          row_indices_->push_back(rr);
          ++index;
        }
      } // rr
      (*column_pointers_)[cc + 1] = index;
    } // cc
  } // CommonSparseMatrix(...)

  template <class OtherMatrixType>
  explicit CommonSparseMatrix(
      const OtherMatrixType& mat,
      const typename std::enable_if<Common::MatrixAbstraction<OtherMatrixType>::is_matrix, bool>::type prune,
      const size_t num_mutexes)
    : CommonSparseMatrix(mat,
                         prune,
                         Common::FloatCmp::DefaultEpsilon<ScalarType>::value()
                             / (1000. * Common::MatrixAbstraction<OtherMatrixType>::cols(mat)),
                         num_mutexes)
  {
  } // CommonSparseMatrix(...)

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      num_rows_ = other.num_rows_;
      num_cols_ = other.num_cols_;
      entries_ = other.unshareable_ ? std::make_shared<EntriesVectorType>(*other.entries_) : other.entries_;
      column_pointers_ = other.column_pointers_;
      row_indices_ = other.row_indices_;
      mutexes_ =
          other.unshareable_ ? std::make_shared<std::vector<std::mutex>>(other.mutexes_->size()) : other.mutexes_;
      eps_ = other.eps_;
      unshareable_ = false;
    }
    return *this;
  }

  template <class OtherMatrixImp>
  typename std::enable_if_t<XT::Common::MatrixAbstraction<OtherMatrixImp>::is_matrix, ThisType>&
  assign(const OtherMatrixImp& other, const SparsityPatternDefault& pattern)
  {

    using MatAbstrType = Common::MatrixAbstraction<OtherMatrixImp>;
    clear();
    num_rows_ = MatAbstrType::rows(other);
    num_cols_ = MatAbstrType::cols(other);
    for (size_t cc = 0; cc < num_cols_; ++cc) {
      for (size_t rr = 0; rr < num_rows_; ++rr) {
        if (XT::Common::FloatCmp::ne(MatAbstrType::get_entry(other, rr, cc), ScalarType(0), 0., eps_ / num_cols_)
            && std::find(pattern.inner(rr).begin(), pattern.inner(rr).end(), cc) != pattern.inner(rr).end()) {
          entries_->push_back(MatAbstrType::get_entry(other, rr, cc));
          row_indices_->push_back(rr);
        }
      } // rr
      (*column_pointers_)[cc + 1] = row_indices_->size();
    } // cc
    return *this;
  }

  template <class OtherMatrixImp>
  typename std::enable_if_t<XT::Common::MatrixAbstraction<OtherMatrixImp>::is_matrix, ThisType>&
  operator=(const OtherMatrixImp& other)
  {
    assign(other, dense_pattern(num_rows_, num_cols_));
  }


  void deep_copy(const ThisType& other)
  {
    ensure_uniqueness();
    num_rows_ = other.num_rows_;
    num_cols_ = other.num_cols_;
    *entries_ = *other.entries_;
    *column_pointers_ = *other.column_pointers_;
    *row_indices_ = *other.row_indices_;
  }

  void clear()
  {
    ensure_uniqueness();
    entries_->clear();
    std::fill(column_pointers_->begin(), column_pointers_->end(), 0);
    row_indices_->clear();
  }

  /// \name Required by ContainerInterface.
  /// \{
  inline ThisType copy() const
  {
    ThisType ret(*this);
    ret.entries_ = std::make_shared<EntriesVectorType>(*entries_);
    ret.column_pointers_ = std::make_shared<IndexVectorType>(*column_pointers_);
    ret.row_indices_ = std::make_shared<IndexVectorType>(*row_indices_);
    ret.mutexes_ = std::make_shared<std::vector<std::mutex>>(mutexes_->size());
    return ret;
  }

  inline void scal(const ScalarType& alpha)
  {
    ensure_uniqueness();
    const internal::VectorLockGuard DUNE_UNUSED(guard)(*mutexes_);
    std::transform(
        entries_->begin(), entries_->end(), entries_->begin(), std::bind1st(std::multiplies<ScalarType>(), alpha));
  }

  inline void axpy(const ScalarType& alpha, const ThisType& xx)
  {
    ensure_uniqueness();
    const internal::VectorLockGuard DUNE_UNUSED(guard)(*mutexes_);
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

  //! Matrix-Vector multiplication for arbitrary vectors that support operator[]
  template <class XX, class YY>
  inline std::enable_if_t<!std::is_base_of<CommonSparseVector<ScalarType>, XX>::value
                              && !std::is_base_of<CommonSparseVector<ScalarType>, YY>::value
                              && XT::Common::VectorAbstraction<XX>::is_vector
                              && XT::Common::VectorAbstraction<YY>::is_vector,
                          void>
  mv(const XX& xx, YY& yy) const
  {
    std::fill(yy.begin(), yy.end(), ScalarType(0));
    const auto& entries = *entries_;
    const auto& column_pointers = *column_pointers_;
    const auto& row_indices = *row_indices_;
    for (size_t cc = 0; cc < num_cols_; ++cc) {
      const size_t end = column_pointers[cc + 1];
      for (size_t kk = column_pointers[cc]; kk < end; ++kk)
        yy[row_indices[kk]] += entries[kk] * xx[cc];
    }
  }

  void mv(const CommonSparseVector<ScalarType>& xx, CommonSparseVector<ScalarType>& yy) const
  {
    yy.clear();
    const auto& entries = *entries_;
    const auto& column_pointers = *column_pointers_;
    const auto& row_indices = *row_indices_;
    const auto& vec_entries = xx.entries();
    const auto& vec_indices = xx.indices();
    thread_local std::vector<ScalarType> tmp_vec;
    tmp_vec.resize(num_rows_);
    std::fill(tmp_vec.begin(), tmp_vec.end(), 0.);
    for (size_t ii = 0; ii < vec_entries.size(); ++ii) {
      const size_t cc = vec_indices[ii];
      const size_t end = column_pointers[cc + 1];
      for (size_t kk = column_pointers[cc]; kk < end; ++kk)
        tmp_vec[row_indices[kk]] += entries[kk] * vec_entries[ii];
    }
    for (size_t cc = 0; cc < num_cols_; ++cc)
      if (XT::Common::FloatCmp::ne(tmp_vec[cc], ScalarType(0.), 0., eps_))
        yy.set_new_entry(cc, tmp_vec[cc]);
  }

  //! TransposedMatrix-Vector multiplication for arbitrary vectors that support operator[]
  template <class XX, class YY>
  inline std::enable_if_t<!std::is_base_of<CommonSparseVector<ScalarType>, XX>::value
                              && !std::is_base_of<CommonSparseVector<ScalarType>, YY>::value
                              && XT::Common::VectorAbstraction<XX>::is_vector
                              && XT::Common::VectorAbstraction<YY>::is_vector,
                          void>
  mtv(const XX& xx, YY& yy) const
  {
    std::fill(yy.begin(), yy.end(), ScalarType(0));
    const auto& entries = *entries_;
    const auto& column_pointers = *column_pointers_;
    const auto& row_indices = *row_indices_;
    for (size_t cc = 0; cc < num_cols_; ++cc) {
      const size_t end = column_pointers[cc + 1];
      for (size_t kk = column_pointers[cc]; kk < end; ++kk)
        yy[cc] += entries[kk] * xx[row_indices[kk]];
    }
  }

  void mtv(const CommonSparseVector<ScalarType>& xx, CommonSparseVector<ScalarType>& yy) const
  {
    yy.clear();
    const auto& entries = *entries_;
    const auto& column_pointers = *column_pointers_;
    const auto& row_indices = *row_indices_;
    for (size_t cc = 0; cc < num_cols_; ++cc) {
      const size_t end = column_pointers[cc + 1];
      ScalarType new_entry(0);
      for (size_t kk = column_pointers[cc]; kk < end; ++kk)
        new_entry += entries[kk] * xx.get_entry(row_indices[kk]);
      if (XT::Common::FloatCmp::ne(new_entry, ScalarType(0.), 0., eps_))
        yy.set_new_entry(cc, new_entry);
    }
  }

  inline void add_to_entry(const size_t rr, const size_t cc, const ScalarType& value)
  {
    ensure_uniqueness();
    internal::LockGuard DUNE_UNUSED(lock)(*mutexes_, rr, rows());
    entries_->operator[](get_entry_index(rr, cc)) += value;
  }

  inline ScalarType get_entry(const size_t rr, const size_t cc) const
  {
    const size_t index = get_entry_index(rr, cc, false);
    return index == size_t(-1) ? ScalarType(0) : (*entries_)[index];
  }

  inline void set_entry(const size_t rr, const size_t cc, const ScalarType value)
  {
    ensure_uniqueness();
    entries_->operator[](get_entry_index(rr, cc)) = value;
  }

  inline void start_column()
  {
    ensure_uniqueness();
    if (column_pointers_->empty())
      column_pointers_->push_back(0);
  }

  inline void end_column()
  {
    ensure_uniqueness();
    column_pointers_->push_back(row_indices_->size());
  }

  inline void push_entry(const size_t cc, const ScalarType value)
  {
    ensure_uniqueness();
    entries_->push_back(value);
    row_indices_->push_back(cc);
  }

  inline void clear_row(const size_t rr)
  {
    ensure_uniqueness();
    for (size_t kk = 0; kk < entries_->size(); ++kk) {
      if ((*row_indices_)[kk] == rr)
        (*entries_)[kk] = ScalarType(0);
    }
  }

  inline void clear_col(const size_t cc)
  {
    ensure_uniqueness();
    std::fill(
        entries_->begin() + (*column_pointers_)[cc], entries_->begin() + (*column_pointers_)[cc + 1], ScalarType(0));
  }

  inline void unit_row(const size_t rr)
  {
    clear_row(rr);
    set_entry(rr, rr, ScalarType(1));
  }

  inline void unit_col(const size_t cc)
  {
    clear_col(cc);
    set_entry(cc, cc, ScalarType(1));
  }

  bool valid() const
  {
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
                                         const EpsType eps = Common::FloatCmp::DefaultEpsilon<ScalarType>::value()
                                                             / 1000.) const override
  {
    SparsityPatternDefault ret(num_rows_);
    for (size_t cc = 0; cc < num_cols_; ++cc) {
      for (size_t kk = (*column_pointers_)[cc]; kk < (*column_pointers_)[cc + 1]; ++kk) {
        if (!prune || !is_zero((*entries_)[kk], eps))
          ret.insert((*row_indices_)[kk], cc);
      }
    } // cc
    return ret;
  } // ... pattern(...)

  /// \}

  template <class OtherMatrixImp>
  typename std::enable_if_t<XT::Common::MatrixAbstraction<OtherMatrixImp>::is_matrix
                                && !(std::is_base_of<ThisType, OtherMatrixImp>::value),
                            void>
  rightmultiply(const OtherMatrixImp& other)
  {
    ensure_uniqueness();
    EntriesVectorType new_entries;
    new_entries.reserve(entries_->size());
    IndexVectorType new_column_pointers(num_cols_ + 1, 0);
    IndexVectorType new_row_indices;
    new_row_indices.reserve(row_indices_->size());
    ScalarType new_entry(0);
    size_t index = 0;
    for (size_t cc = 0; cc < num_cols_; ++cc) {
      for (size_t rr = 0; rr < num_rows_; ++rr) {
        new_entry = 0;
        for (size_t col = 0; col < num_cols_; ++col) {
          for (size_t kk = (*column_pointers_)[col]; kk < (*column_pointers_)[col + 1]; ++kk)
            if ((*row_indices_)[kk] == rr)
              new_entry += (*entries_)[kk] * XT::Common::MatrixAbstraction<OtherMatrixImp>::get_entry(other, col, cc);
        } // col
        if (XT::Common::FloatCmp::ne(new_entry, 0., 0., eps_ / num_cols_)) {
          new_entries.push_back(new_entry);
          new_row_indices.push_back(rr);
          ++index;
        }
      } // rr
      new_column_pointers[cc + 1] = index;
    } // cc
    *entries_ = new_entries;
    *column_pointers_ = new_column_pointers;
    *row_indices_ = new_row_indices;
  } // void rightmultiply(...)

  void rightmultiply(const ThisType& other)
  {
    ensure_uniqueness();
    thread_local EntriesVectorType new_entries;
    thread_local IndexVectorType new_column_pointers(num_cols_ + 1, 0);
    thread_local IndexVectorType new_row_indices;
    new_entries.clear();
    new_row_indices.clear();
    const auto& other_entries = *other.entries_;
    const auto& other_column_pointers = *other.column_pointers_;
    const auto& other_row_indices = *other.row_indices_;
    thread_local std::vector<ScalarType> dense_column(num_rows_, 0.);
    dense_column.resize(num_rows_);
    for (size_t cc = 0; cc < other.cols(); ++cc) {
      std::fill(dense_column.begin(), dense_column.end(), 0.);
      size_t other_col_end = other_column_pointers[cc + 1];
      for (size_t kk = other_column_pointers[cc]; kk < other_col_end; ++kk) {
        size_t col = other_row_indices[kk];
        size_t col_end = (*column_pointers_)[col + 1];
        for (size_t ll = (*column_pointers_)[col]; ll < col_end; ++ll)
          dense_column[(*row_indices_)[ll]] += (*entries_)[ll] * other_entries[kk];
      } // kk
      for (size_t rr = 0; rr < num_rows_; ++rr) {
        if (XT::Common::FloatCmp::ne(dense_column[rr], ScalarType(0.), 0., eps_ / num_cols_)) {
          new_entries.push_back(dense_column[rr]);
          new_row_indices.push_back(rr);
        }
      } // ii
      new_column_pointers[cc + 1] = new_row_indices.size();
    } // cc
    *entries_ = new_entries;
    *column_pointers_ = new_column_pointers;
    *row_indices_ = new_row_indices;
  } // void rightmultiply(...)

  using InterfaceType::operator+;
  using InterfaceType::operator-;
  using InterfaceType::operator+=;
  using InterfaceType::operator-=;

  ScalarType* entries()
  {
    ensure_uniqueness();
    unshareable_ = true;
    return entries_->data();
  }

  const ScalarType* entries() const
  {
    ensure_uniqueness();
    unshareable_ = true;
    return entries_->data();
  }

  size_t* outer_index_ptr()
  {
    ensure_uniqueness();
    unshareable_ = true;
    return column_pointers_->data();
  }

  const size_t* outer_index_ptr() const
  {
    ensure_uniqueness();
    unshareable_ = true;
    return column_pointers_->data();
  }

  size_t* inner_index_ptr()
  {
    ensure_uniqueness();
    unshareable_ = true;
    return row_indices_->data();
  }

  const size_t* inner_index_ptr() const
  {
    ensure_uniqueness();
    unshareable_ = true;
    return row_indices_->data();
  }

  const EpsType& eps()
  {
    return eps_;
  }

  inline void ensure_uniqueness() const
  {
    if (!entries_.unique()) {
      assert(!unshareable_);
      const internal::VectorLockGuard DUNE_UNUSED(guard)(*mutexes_);
      if (!entries_.unique()) {
        entries_ = std::make_shared<EntriesVectorType>(*entries_);
        column_pointers_ = std::make_shared<IndexVectorType>(*column_pointers_);
        row_indices_ = std::make_shared<IndexVectorType>(*row_indices_);
        mutexes_ = std::make_shared<std::vector<std::mutex>>(mutexes_->size());
      }
    }
  } // ... ensure_uniqueness(...)

private:
  size_t get_entry_index(const size_t rr, const size_t cc, const bool throw_if_not_in_pattern = true) const
  {
    const auto& column_offset = column_pointers_->operator[](cc);
    const auto& next_column_offset = column_pointers_->operator[](cc + 1);
    const auto it = row_indices_->begin() + column_offset;
    const auto it_end = row_indices_->begin() + next_column_offset;
    const auto entry_it = std::lower_bound(it, it_end, rr);
    if (entry_it != it_end && *entry_it == rr)
      return column_offset + std::distance(it, entry_it);
    if (throw_if_not_in_pattern)
      DUNE_THROW(Common::Exceptions::index_out_of_range, "Entry is not in the sparsity pattern!");
    return size_t(-1);
  }

  bool is_zero(const ScalarType val, const ScalarType eps) const
  {
    const auto factor = (num_cols_ == 0 ? 1. : ScalarType(num_cols_));
    const ScalarType tol = eps / factor;
    return XT::Common::FloatCmp::eq(val, ScalarType(0.), 0., tol);
  }

  size_t num_rows_, num_cols_;
  mutable std::shared_ptr<EntriesVectorType> entries_;
  mutable std::shared_ptr<IndexVectorType> column_pointers_;
  mutable std::shared_ptr<IndexVectorType> row_indices_;
  mutable std::shared_ptr<std::vector<std::mutex>> mutexes_;
  mutable EpsType eps_;
  mutable bool unshareable_;
}; // class CommonSparseMatrix<..., Common::StorageLayout::csc>

/**
 * \brief A matrix implementation checking whether the matrix is sparse enough to use sparse matrix operations.
 */
template <class DenseMatrixImp, class SparseMatrixImp>
class CommonSparseOrDenseMatrix
    : public MatrixInterface<internal::CommonSparseOrDenseMatrixTraits<DenseMatrixImp, SparseMatrixImp>,
                             typename SparseMatrixImp::ScalarType>
{
  using ThisType = CommonSparseOrDenseMatrix;
  using InterfaceType = MatrixInterface<internal::CommonSparseOrDenseMatrixTraits<DenseMatrixImp, SparseMatrixImp>,
                                        typename SparseMatrixImp::ScalarType>;

public:
  using typename InterfaceType::RealType;
  using typename InterfaceType::ScalarType;
  using typename InterfaceType::Traits;
  using EpsType = typename Traits::EpsType;
  using DenseMatrixType = typename Traits::DenseMatrixType;
  using SparseMatrixType = typename Traits::SparseMatrixType;
  static constexpr double sparse_limit = 0.1;

  /**
  * \brief This is the constructor of interest which creates a sparse matrix.
  */
  CommonSparseOrDenseMatrix(const size_t rr,
                            const size_t cc,
                            const SparsityPatternDefault& patt,
                            const size_t num_mutexes = 1,
                            EpsType eps = Common::FloatCmp::DefaultEpsilon<ScalarType>::value() / 1000.)
    : num_rows_(rr)
    , num_cols_(cc)
  {
    size_t nnz = 0;
    for (size_t row = 0; row < num_rows_; ++row)
      nnz += patt.inner(row).size();
    size_t num_entries = rr * cc;
    const double density = double(nnz) / double(num_entries);
    sparse_ = density < sparse_limit;
    if (sparse_) {
      sparse_matrix_ = SparseMatrixType(rr, cc, patt, num_mutexes, eps);
      dense_matrix_ = DenseMatrixType(0, 0, patt, num_mutexes);
    } else {
      sparse_matrix_ = SparseMatrixType(0, 0, patt, num_mutexes, eps);
      dense_matrix_ = DenseMatrixType(rr, cc, patt, num_mutexes);
    }
  } // CommonSparseOrDenseMatrix(rr, cc, patt, num_mutexes)

  CommonSparseOrDenseMatrix(const size_t rr = 0,
                            const size_t cc = 0,
                            const ScalarType& value = ScalarType(0),
                            const size_t num_mutexes = 1,
                            bool use_sparse_if_zero = true,
                            const EpsType eps = Common::FloatCmp::DefaultEpsilon<ScalarType>::value() / 1000.)
    : num_rows_(rr)
    , num_cols_(cc)
  {
    if (!is_zero(value, eps) || !use_sparse_if_zero) {
      sparse_matrix_ = SparseMatrixType(0, 0, value, num_mutexes, eps);
      dense_matrix_ = DenseMatrixType(rr, cc, value, num_mutexes);
    } else {
      sparse_matrix_ = SparseMatrixType(rr, cc, value, num_mutexes, eps);
      dense_matrix_ = DenseMatrixType(0, 0, value, num_mutexes);
      sparse_ = true;
    }
  }

  CommonSparseOrDenseMatrix(const size_t rr, const size_t cc, const size_t num_mutexes, bool use_sparse = true)
    : CommonSparseOrDenseMatrix(rr, cc, ScalarType(0.), num_mutexes, use_sparse)
  {
  }

  CommonSparseOrDenseMatrix(const ThisType& other)
    : num_rows_(other.num_rows_)
    , num_cols_(other.num_cols_)
    , sparse_(other.sparse_)
    , sparse_matrix_(other.sparse_matrix_)
    , dense_matrix_(other.dense_matrix_)
  {
  }

  template <class OtherMatrixType>
  explicit CommonSparseOrDenseMatrix(
      const OtherMatrixType& mat,
      const typename std::enable_if<Common::MatrixAbstraction<OtherMatrixType>::is_matrix, bool>::type prune = true,
      const EpsType eps_in = Common::FloatCmp::DefaultEpsilon<ScalarType>::value() / 1000.,
      const size_t num_mutexes = 1)
    : num_rows_(Common::MatrixAbstraction<OtherMatrixType>::rows(mat))
    , num_cols_(Common::MatrixAbstraction<OtherMatrixType>::cols(mat))
  {
    auto eps = eps_in / Common::MatrixAbstraction<OtherMatrixType>::cols(mat);
    // check sparsity (if prune = false, the sparsity may not be checked correctly)
    size_t nnz = 0.;
    for (size_t rr = 0; rr < num_rows_; ++rr)
      for (size_t cc = 0; cc < num_cols_; ++cc) {
        const auto& value = Common::MatrixAbstraction<OtherMatrixType>::get_entry(mat, rr, cc);
        nnz += XT::Common::FloatCmp::ne(value, ScalarType(0), 0., eps);
      }
    const double density = double(nnz) / double(num_rows_ * num_cols_);
    sparse_ = density < sparse_limit;
    if (sparse_) {
      sparse_matrix_ = SparseMatrixType(mat, prune, eps_in, num_mutexes);
      dense_matrix_ = DenseMatrixType(0, 0, num_mutexes);
    } else {
      sparse_matrix_ = SparseMatrixType(0, 0, num_mutexes);
      dense_matrix_ = DenseMatrixType(mat, num_mutexes);
    } // else (sparse_)
  } // CommonSparseOrDenseMatrix(...)

  template <class OtherMatrixType>
  explicit CommonSparseOrDenseMatrix(
      const OtherMatrixType& mat,
      const typename std::enable_if<Common::MatrixAbstraction<OtherMatrixType>::is_matrix, bool>::type prune,
      const size_t num_mutexes)
    : CommonSparseOrDenseMatrix(mat,
                                prune,
                                Common::FloatCmp::DefaultEpsilon<ScalarType>::value()
                                    / (1000. * Common::MatrixAbstraction<OtherMatrixType>::cols(mat)),
                                num_mutexes)
  {
  } // CommonSparseOrDenseMatrix(...)

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      num_rows_ = other.num_rows_;
      num_cols_ = other.num_cols_;
      sparse_ = other.sparse_;
      sparse_matrix_ = other.sparse_matrix_;
      dense_matrix_ = other.dense_matrix_;
    }
    return *this;
  }

  ThisType& operator=(const DenseMatrixType& other)
  {
    size_t nnz = 0.;
    auto tol = sparse_matrix_.eps() / Common::MatrixAbstraction<DenseMatrixType>::cols(other);
    for (size_t rr = 0; rr < num_rows_; ++rr)
      for (size_t cc = 0; cc < num_cols_; ++cc) {
        const auto& value = Common::MatrixAbstraction<DenseMatrixType>::get_entry(other, rr, cc);
        nnz += XT::Common::FloatCmp::ne(value, ScalarType(0), 0., tol);
      }
    const double density = double(nnz) / (Common::MatrixAbstraction<DenseMatrixType>::rows(other)
                                          * Common::MatrixAbstraction<DenseMatrixType>::cols(other));
    sparse_ = density < sparse_limit;
    if (sparse_)
      sparse_matrix_ = other;
    else
      dense_matrix_ = other;
    return *this;
  }

  void deep_copy(const ThisType& other)
  {
    num_rows_ = other.num_rows_;
    num_cols_ = other.num_cols_;
    sparse_ = other.sparse_;
    sparse_matrix_.deep_copy(other.sparse_matrix_);
    dense_matrix_.deep_copy(other.dense_matrix_);
  }

  /// \name Required by ContainerInterface.
  /// \{
  inline ThisType copy() const
  {
    ThisType ret(*this);
    ret.sparse_matrix_ = sparse_matrix_.copy();
    ret.dense_matrix_ = dense_matrix_.copy();
    return ret;
  }

  inline void scal(const ScalarType& alpha)
  {
    sparse_ ? sparse_matrix_.scal(alpha) : dense_matrix_.scal(alpha);
  }

  inline void axpy(const ScalarType& alpha, const ThisType& xx)
  {
    assert(sparse_ == xx.sparse_);
    sparse_ ? sparse_matrix_.axpy(alpha, xx.sparse_matrix_) : dense_matrix_.axpy(alpha, xx.dense_matrix_);
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

  //! Matrix-Vector multiplication for arbitrary vectors that support operator[]
  template <class XX, class YY>
  inline std::enable_if_t<XT::Common::VectorAbstraction<XX>::is_vector && XT::Common::VectorAbstraction<YY>::is_vector,
                          void>
  mv(const XX& xx, YY& yy) const
  {
    sparse_ ? sparse_matrix_.mv(xx, yy) : dense_matrix_.mv(xx, yy);
  }

  //! TransposedMatrix-Vector multiplication for arbitrary vectors that support operator[]
  template <class XX, class YY>
  inline std::enable_if_t<XT::Common::VectorAbstraction<XX>::is_vector && XT::Common::VectorAbstraction<YY>::is_vector,
                          void>
  mtv(const XX& xx, YY& yy) const
  {
    sparse_ ? sparse_matrix_.mtv(xx, yy) : dense_matrix_.mtv(xx, yy);
  }

  inline void add_to_entry(const size_t rr, const size_t cc, const ScalarType& value)
  {
    sparse_ ? sparse_matrix_.add_to_entry(rr, cc, value) : dense_matrix_.add_to_entry(rr, cc, value);
  }

  inline ScalarType get_entry(const size_t rr, const size_t cc) const
  {
    return sparse_ ? sparse_matrix_.get_entry(rr, cc) : dense_matrix_.get_entry(rr, cc);
  }

  inline void set_entry(const size_t rr, const size_t cc, const ScalarType value)
  {
    sparse_ ? sparse_matrix_.set_entry(rr, cc, value) : dense_matrix_.set_entry(rr, cc, value);
  }

  inline void clear_row(const size_t rr)
  {
    sparse_ ? sparse_matrix_.clear_row(rr) : dense_matrix_.clear_row(rr);
  }

  inline void clear_col(const size_t cc)
  {
    sparse_ ? sparse_matrix_.clear_col(cc) : dense_matrix_.clear_col(cc);
  }

  inline void unit_row(const size_t rr)
  {
    clear_row(rr);
    set_entry(rr, rr, ScalarType(1));
  }

  inline void unit_col(const size_t cc)
  {
    clear_col(cc);
    set_entry(cc, cc, ScalarType(1));
  }

  bool valid() const
  {
    return sparse_ ? sparse_matrix_.valid() : dense_matrix_.valid();
  }

  virtual size_t non_zeros() const override final
  {
    return sparse_ ? sparse_matrix_.non_zeros() : dense_matrix_.non_zeros();
  }

  virtual SparsityPatternDefault
  pattern(const bool prune = false,
          const EpsType eps = Common::FloatCmp::DefaultEpsilon<ScalarType>::value()) const override
  {
    return sparse_ ? sparse_matrix_.pattern(prune, eps) : dense_matrix_.pattern(prune, eps);
  } // ... pattern(...)

  /// \}

  template <class DuneDenseMatrixImp>
  void copy_to_densematrix(DuneDenseMatrixImp& ret) const
  {
    sparse_ ? sparse_matrix_.copy_to_densematrix(ret) : dense_matrix_.copy_to_densematrix(ret);
  }

  bool sparse() const
  {
    return sparse_;
  }

  SparseMatrixType& sparse_matrix()
  {
    return sparse_matrix_;
  }

  const SparseMatrixType& sparse_matrix() const
  {
    return sparse_matrix_;
  }

  DenseMatrixType& dense_matrix()
  {
    return dense_matrix_;
  }

  const DenseMatrixType& dense_matrix() const
  {
    return dense_matrix_;
  }

  template <class MatrixType>
  void rightmultiply(const MatrixType& other)
  {
    sparse_ ? sparse_matrix_.rightmultiply(other) : dense_matrix_.rightmultiply(other);
  } // void rightmultiply(...)

  void rightmultiply(const ThisType& other)
  {
    if (other.sparse())
      sparse_ ? sparse_matrix_.rightmultiply(other.sparse_matrix())
              : dense_matrix_.rightmultiply(other.sparse_matrix());
    else
      sparse_ ? sparse_matrix_.rightmultiply(other.dense_matrix()) : dense_matrix_.rightmultiply(other.dense_matrix());
  } // void rightmultiply(...)

  using InterfaceType::operator+;
  using InterfaceType::operator-;
  using InterfaceType::operator+=;
  using InterfaceType::operator-=;

  // clang does not find the conversion operators from the interface if they are not redefined here
  template <int ROWS, int COLS>
  explicit operator Dune::FieldMatrix<ScalarType, ROWS, COLS>() const
  {
    return InterfaceType::operator Dune::FieldMatrix<ScalarType, ROWS, COLS>();
  }

  template <int ROWS, int COLS>
  explicit operator std::unique_ptr<Dune::FieldMatrix<ScalarType, ROWS, COLS>>() const
  {
    return InterfaceType::operator std::unique_ptr<Dune::FieldMatrix<ScalarType, ROWS, COLS>>();
  }

  template <int ROWS, int COLS>
  explicit operator std::unique_ptr<XT::Common::FieldMatrix<ScalarType, ROWS, COLS>>() const
  {
    return InterfaceType::operator std::unique_ptr<XT::Common::FieldMatrix<ScalarType, ROWS, COLS>>();
  }

  explicit operator Dune::DynamicMatrix<ScalarType>() const
  {
    return InterfaceType::operator Dune::DynamicMatrix<ScalarType>();
  }

private:
  bool is_zero(const ScalarType val, const ScalarType eps) const
  {
    const auto factor = (num_cols_ == 0 ? 1. : ScalarType(num_cols_));
    const ScalarType tol = eps / factor;
    return XT::Common::FloatCmp::eq(val, ScalarType(0.), 0., tol);
  }

  size_t num_rows_, num_cols_;
  bool sparse_;
  mutable SparseMatrixType sparse_matrix_;
  mutable DenseMatrixType dense_matrix_;
}; // class CommonSparseOrDenseMatrix<...>

template <class ScalarType = double>
using CommonSparseMatrixCsr = CommonSparseMatrix<ScalarType, Common::StorageLayout::csr>;

template <class ScalarType = double>
using CommonSparseMatrixCsc = CommonSparseMatrix<ScalarType, Common::StorageLayout::csc>;

template <class ScalarType = double>
using CommonSparseOrDenseMatrixCsr =
    CommonSparseOrDenseMatrix<CommonDenseMatrix<ScalarType>, CommonSparseMatrixCsr<ScalarType>>;

template <class ScalarType = double>
using CommonSparseOrDenseMatrixCsc =
    CommonSparseOrDenseMatrix<CommonDenseMatrix<ScalarType>, CommonSparseMatrixCsc<ScalarType>>;


} // namespace LA
namespace Common {


template <class T>
struct MatrixAbstraction<LA::CommonSparseMatrixCsr<T>>
    : public LA::internal::MatrixAbstractionBase<LA::CommonSparseMatrixCsr<T>>
{
  using BaseType = LA::internal::MatrixAbstractionBase<LA::CommonSparseMatrixCsr<T>>;

  template <size_t rows = BaseType::static_rows, size_t cols = BaseType::static_cols, class FieldType = T>
  using MatrixTypeTemplate = LA::CommonSparseMatrixCsr<FieldType>;

  static const constexpr Common::StorageLayout storage_layout = Common::StorageLayout::csr;
};

template <class T>
struct MatrixAbstraction<LA::CommonSparseMatrixCsc<T>>
    : public LA::internal::MatrixAbstractionBase<LA::CommonSparseMatrixCsc<T>>
{
  using BaseType = LA::internal::MatrixAbstractionBase<LA::CommonSparseMatrixCsc<T>>;

  template <size_t rows = BaseType::static_rows, size_t cols = BaseType::static_cols, class FieldType = T>
  using MatrixTypeTemplate = LA::CommonSparseMatrixCsc<FieldType>;

  static const constexpr Common::StorageLayout storage_layout = Common::StorageLayout::csc;
};

template <class DenseMatrixImp, class SparseMatrixImp>
struct MatrixAbstraction<LA::CommonSparseOrDenseMatrix<DenseMatrixImp, SparseMatrixImp>>
    : public LA::internal::MatrixAbstractionBase<LA::CommonSparseOrDenseMatrix<DenseMatrixImp, SparseMatrixImp>>
{
  using BaseType = LA::internal::MatrixAbstractionBase<LA::CommonSparseOrDenseMatrix<DenseMatrixImp, SparseMatrixImp>>;

  template <size_t rows = BaseType::static_rows,
            size_t cols = BaseType::static_cols,
            class FieldType = typename BaseType::ScalarType>
  using MatrixTypeTemplate = LA::CommonSparseOrDenseMatrix<
      typename MatrixAbstraction<DenseMatrixImp>::template MatrixTypeTemplate<rows, cols, FieldType>,
      typename MatrixAbstraction<SparseMatrixImp>::template MatrixTypeTemplate<rows, cols, FieldType>>;

  static const constexpr Common::StorageLayout storage_layout = Common::StorageLayout::other;
};


} // namespace Common
} // namespace XT
} // namespace Dune


// begin: this is what we need for the lib
#if DUNE_XT_WITH_PYTHON_BINDINGS


extern template class Dune::XT::LA::CommonSparseMatrix<double, Dune::XT::Common::StorageLayout::csr>;
extern template class Dune::XT::LA::CommonSparseMatrix<double, Dune::XT::Common::StorageLayout::csc>;


#endif // DUNE_XT_WITH_PYTHON_BINDINGS
// end: this is what we need for the lib


#endif // DUNE_XT_LA_CONTAINER_COMMON_MATRIX_SPARSE_HH
