// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// The copyright lies with the authors of this file (see below).
// License: Dual licensed as  BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Barbara Verfürth (2015)
//   Felix Schindler  (2013 - 2016)
//   Rene Milk        (2014 - 2015)
//   Tobias Leibner   (2014)

#ifndef DUNE_XT_LA_CONTAINER_COMMON_HH
#define DUNE_XT_LA_CONTAINER_COMMON_HH

#include <cmath>
#include <initializer_list>
#include <memory>
#include <type_traits>
#include <vector>
#include <complex>

#include <dune/common/dynvector.hh>
#include <dune/common/dynmatrix.hh>
#include <dune/common/densematrix.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/ftraits.hh>

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

namespace internal {

/// Traits for CommonDenseVector
template <class ScalarImp = double>
class CommonDenseVectorTraits
{
public:
  typedef typename Dune::FieldTraits<ScalarImp>::field_type ScalarType;
  typedef typename Dune::FieldTraits<ScalarImp>::real_type RealType;
  typedef CommonDenseVector<ScalarType> derived_type;
  typedef Dune::DynamicVector<ScalarType> BackendType;
};

template <class ScalarImp = double>
class CommonDenseMatrixTraits
{
public:
  typedef typename Dune::FieldTraits<ScalarImp>::field_type ScalarType;
  typedef typename Dune::FieldTraits<ScalarImp>::real_type RealType;
  typedef CommonDenseMatrix<ScalarType> derived_type;
  typedef Dune::DynamicMatrix<ScalarType> BackendType;
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

  CommonDenseVector(const ThisType& other) = default;

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
    backend_ = std::make_shared<BackendType>(other);
    return *this;
  }

  /// \name Required by the ProvidesBackend interface.
  /// \{

  BackendType& backend()
  {
    ensure_uniqueness();
    return *backend_;
  } // ... backend(...)

  const BackendType& backend() const
  {
    ensure_uniqueness();
    return *backend_;
  } // ... backend(...)

  /// \}
  /// \name Required by ProvidesDataAccess.
  /// \{

  ScalarType* data()
  {
    return &(backend()[0]);
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
    backend() *= alpha;
  } // ... scal(...)

  void axpy(const ScalarType& alpha, const ThisType& xx)
  {
    if (xx.size() != size())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of x (" << xx.size() << ") does not match the size of this (" << size() << ")!");
    ensure_uniqueness();
    auto& this_ref = *backend_;
    const auto& xx_ref = *(xx.backend_);
    for (size_t ii = 0; ii < this_ref.size(); ++ii)
      this_ref[ii] += alpha * xx_ref[ii];
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
    assert(ii < size());
    ensure_uniqueness();
    backend_->operator[](ii) += value;
  } // ... add_to_entry(...)

  void set_entry(const size_t ii, const ScalarType& value)
  {
    assert(ii < size());
    ensure_uniqueness();
    backend_->operator[](ii) = value;
  } // ... set_entry(...)

  ScalarType get_entry(const size_t ii) const
  {
    assert(ii < size());
    return backend_->operator[](ii);
  } // ... get_entry(...)

private:
  inline ScalarType& get_entry_ref(const size_t ii)
  {
    return backend()[ii];
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
    if (other.size() != size())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    return backend_->operator*(*(other.backend_));
  } // ... dot(...)

  virtual RealType l1_norm() const override final
  {
    return backend_->one_norm();
  }

  virtual RealType l2_norm() const override final
  {
    return backend_->two_norm();
  }

  virtual RealType sup_norm() const override final
  {
    return backend_->infinity_norm();
  }

  virtual void add(const ThisType& other, ThisType& result) const override final
  {
    if (other.size() != size())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    if (result.size() != size())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of result (" << result.size() << ") does not match the size of this (" << size() << ")!");
    BackendType& result_ref = result.backend();
    for (size_t ii = 0; ii < size(); ++ii)
      result_ref[ii] = backend_->operator[](ii) + other.backend_->operator[](ii);
  } // ... add(...)

  virtual void iadd(const ThisType& other) override final
  {
    if (other.size() != size())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    backend() += *(other.backend_);
  } // ... iadd(...)

  virtual void sub(const ThisType& other, ThisType& result) const override final
  {
    if (other.size() != size())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    if (result.size() != size())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of result (" << result.size() << ") does not match the size of this (" << size() << ")!");
    BackendType& result_ref = result.backend();
    for (size_t ii = 0; ii < size(); ++ii)
      result_ref[ii] = backend_->operator[](ii) - other.backend_->operator[](ii);
  } // ... sub(...)

  virtual void isub(const ThisType& other) override final
  {
    if (other.size() != size())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    backend() -= *(other.backend_);
  } // ... isub(...)

  /// \}
  /// \name Imported from VectorInterface.
  /// \{

  using VectorInterfaceType::add;
  using VectorInterfaceType::sub;

  /// \}

private:
  /**
   * \see ContainerInterface
   */
  inline void ensure_uniqueness() const
  {
    if (!backend_.unique())
      backend_ = std::make_shared<BackendType>(*backend_);
  } // ... ensure_uniqueness(...)

  friend class VectorInterface<internal::CommonDenseVectorTraits<ScalarType>, ScalarType>;
  friend class CommonDenseMatrix<ScalarType>;

  mutable std::shared_ptr<BackendType> backend_;
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
    backend_ = other.backend_;
    return *this;
  }

  /**
   *  \note Does a deep copy.
   */
  ThisType& operator=(const BackendType& other)
  {
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
    ensure_uniqueness();
    return *backend_;
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
    backend() *= alpha;
  }

  void axpy(const ScalarType& alpha, const ThisType& xx)
  {
    if (!has_equal_shape(xx))
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The shape of xx (" << xx.rows() << "x" << xx.cols() << ") does not match the shape of this ("
                                     << rows()
                                     << "x"
                                     << cols()
                                     << ")!");
    backend().axpy(alpha, *(xx.backend_));
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
    return backend_->rows();
  }

  inline size_t cols() const
  {
    return backend_->cols();
  }

  inline void mv(const VectorInterface<internal::CommonDenseVectorTraits<ScalarType>, ScalarType>& xx,
                 VectorInterface<internal::CommonDenseVectorTraits<ScalarType>, ScalarType>& yy) const
  {
    mv(xx.as_imp(), yy.as_imp());
  }

  inline void mv(const CommonDenseVector<ScalarType>& xx, CommonDenseVector<ScalarType>& yy) const
  {
    backend_->mv(*(xx.backend_), yy.backend());
  }

  void add_to_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    assert(ii < rows());
    assert(jj < cols());
    backend()[ii][jj] += value;
  } // ... add_to_entry(...)

  void set_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    assert(ii < rows());
    assert(jj < cols());
    backend()[ii][jj] = value;
  } // ... set_entry(...)

  ScalarType get_entry(const size_t ii, const size_t jj) const
  {
    assert(ii < rows());
    assert(jj < cols());
    return backend_->operator[](ii)[jj];
  } // ... get_entry(...)

  void clear_row(const size_t ii)
  {
    if (ii >= rows())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given ii (" << ii << ") is larger than the rows of this (" << rows() << ")!");
    backend()[ii] *= ScalarType(0);
  } // ... clear_row(...)

  void clear_col(const size_t jj)
  {
    if (jj >= cols())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given jj (" << jj << ") is larger than the cols of this (" << cols() << ")!");
    BackendType& backend_ref = backend();
    for (size_t ii = 0; ii < rows(); ++ii)
      backend_ref[ii][jj] = ScalarType(0);
  } // ... clear_col(...)

  void unit_row(const size_t ii)
  {
    if (ii >= cols())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given ii (" << ii << ") is larger than the cols of this (" << cols() << ")!");
    if (ii >= rows())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given ii (" << ii << ") is larger than the rows of this (" << rows() << ")!");
    auto& row = backend()[ii];
    for (size_t jj = 0; jj < cols(); ++jj)
      row[jj] = ScalarType(0);
    row[ii] = ScalarType(1);
  } // ... unit_row(...)

  void unit_col(const size_t jj)
  {
    if (jj >= cols())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given jj (" << jj << ") is larger than the cols of this (" << cols() << ")!");
    if (jj >= rows())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given jj (" << jj << ") is larger than the rows of this (" << rows() << ")!");
    ensure_uniqueness();
    for (size_t ii = 0; ii < rows(); ++ii)
      backend_->operator[](ii)[jj] = ScalarType(0);
    backend_->operator[](jj)[jj] = ScalarType(1);
  } // ... unit_col(...)

  bool valid() const
  {
    for (size_t ii = 0; ii < rows(); ++ii) {
      const auto& row_vec = backend_->operator[](ii);
      for (size_t jj = 0; jj < cols(); ++jj) {
        const auto& entry = row_vec[jj];
        if (Common::isnan(entry) || Common::isinf(entry))
          return false;
      }
    }
    return true;
  } // ... valid(...)

  /// \}

private:
  /**
   * \see ContainerInterface
   */
  inline void ensure_uniqueness() const
  {
    if (!backend_.unique())
      backend_ = std::make_shared<BackendType>(*backend_);
  } // ... ensure_uniqueness(...)

  mutable std::shared_ptr<BackendType> backend_;
}; // class CommonDenseMatrix

/**
 * \brief A sparse matrix implementation of the MatrixInterface with row major memory layout.
 */
template <class ScalarImp = double>
class CommonSparseMatrix : public MatrixInterface<internal::CommonSparseMatrixTraits<ScalarImp>, ScalarImp>,
                           public ProvidesBackend<internal::CommonSparseMatrixTraits<ScalarImp>>
{
  typedef CommonSparseMatrix<ScalarImp> ThisType;
  typedef MatrixInterface<internal::CommonSparseMatrixTraits<ScalarImp>, ScalarImp> MatrixInterfaceType;
  static_assert(!std::is_same<DUNE_STUFF_SSIZE_T, int>::value,
                "You have to manually disable the constructor below which uses DUNE_STUFF_SSIZE_T!");

public:
  typedef internal::CommonSparseMatrixTraits<ScalarImp> Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::ScalarType ScalarType;
  typedef typename Traits::RealType RealType;
  typedef std::vector<size_t> IndexVectorType;


  /**
  * \brief This is the constructor of interest which creates a sparse matrix.
  */
  CommonSparseMatrix(const size_t rr, const size_t cc, const SparsityPatternDefault& patt)
    : num_rows_(rr)
    , num_cols_(cc)
    , backend_(std::make_shared<BackendType>())
    , row_pointers_(std::make_shared<IndexVectorType>(num_rows_ + 1))
    , column_indices_(std::make_shared<IndexVectorType>())
  {
    if (num_rows_ > 0 && num_cols_ > 0) {
      if (size_t(patt.size()) != num_rows_)
        DUNE_THROW(Exceptions::shapes_do_not_match,
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
            DUNE_THROW(Exceptions::shapes_do_not_match,
                       "The size of row " << row << " of the pattern does not match the number of columns of this ("
                                          << num_cols_
                                          << ")!");
#endif // NDEBUG
          column_indices_->push_back(columns[col]);
        }
        backend_->resize(column_indices_->size());
      }
    }
  }

  CommonSparseMatrix(const size_t rr = 0, const size_t cc = 0, const ScalarType& value = ScalarType(0))
    : num_rows_(rr)
    , num_cols_(cc)
    , backend_(std::make_shared<BackendType>(num_rows_ * num_cols_, value))
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

  /// This constructor is needed for the python bindings.
  explicit CommonSparseMatrix(const DUNE_STUFF_SSIZE_T rr, const DUNE_STUFF_SSIZE_T cc = 0)
    : CommonSparseMatrix(size_t(rr), size_t(cc), ScalarType(0))
  {
  }

  explicit CommonSparseMatrix(const int rr, const int cc = 0)
    : CommonSparseMatrix(size_t(rr), size_t(cc), ScalarType(0))
  {
    assert(rr >= 0 && cc >= 0);
  }

  CommonSparseMatrix(const ThisType& other)
    : num_rows_(other.num_rows_)
    , num_cols_(other.num_cols_)
    , backend_(other.backend_)
    , row_pointers_(other.row_pointers_)
    , column_indices_(other.column_indices_)
  {
  }

  template <class OtherMatrixType>
  explicit CommonSparseMatrix(
      const OtherMatrixType& mat,
      const typename std::enable_if<DSC::MatrixAbstraction<OtherMatrixType>::is_matrix, bool>::type prune = false,
      const typename Common::FloatCmp::DefaultEpsilon<ScalarType>::Type eps =
          Common::FloatCmp::DefaultEpsilon<ScalarType>::value())
    : num_rows_(DSC::MatrixAbstraction<OtherMatrixType>::rows(mat))
    , num_cols_(DSC::MatrixAbstraction<OtherMatrixType>::cols(mat))
    , backend_(std::make_shared<BackendType>())
    , row_pointers_(std::make_shared<IndexVectorType>(num_rows_ + 1))
    , column_indices_(std::make_shared<IndexVectorType>())
  {
    for (size_t rr = 0; rr < num_rows_; ++rr) {
      size_t num_nonzero_entries_in_row = 0;
      for (size_t cc = 0; cc < num_cols_; ++cc) {
        const auto& value = DSC::MatrixAbstraction<OtherMatrixType>::get_entry(mat, rr, cc);
        if (!prune || DSC::FloatCmp::ne(value, ScalarType(0), eps)) {
          ++num_nonzero_entries_in_row;
          backend_->push_back(value);
          column_indices_->push_back(cc);
        }
      }
      row_pointers_->operator[](rr + 1) = row_pointers_->operator[](rr) + num_nonzero_entries_in_row;
    }
  } // CommonSparseMatrix(...)

  template <int ROWS, int COLS>
  explicit operator Dune::FieldMatrix<ScalarType, ROWS, COLS>() const
  {
    assert(ROWS == num_rows_ && COLS == num_cols_);
    Dune::FieldMatrix<ScalarType, ROWS, COLS> ret(ScalarType(0));
    for (size_t rr = 0; rr < ROWS; ++rr)
      for (size_t kk = row_pointers_->operator[](rr); kk < row_pointers_->operator[](rr + 1); ++kk)
        ret[rr][column_indices_->operator[](kk)] = backend_->operator[](kk);
    return ret;
  }

  ThisType& operator=(const ThisType& other)
  {
    backend_ = other.backend_;
    row_pointers_ = other.row_pointers_;
    column_indices_ = other.column_indices_;
    num_rows_ = other.num_rows_;
    num_cols_ = other.num_cols_;
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

  inline ThisType copy() const
  {
    CommonSparseMatrix ret(*this);
    ensure_uniqueness();
    return ret;
  }

  inline void scal(const ScalarType& alpha)
  {
    ensure_uniqueness();
    std::transform(
        backend_->begin(), backend_->end(), backend_->begin(), std::bind1st(std::multiplies<ScalarType>(), alpha));
  }

  inline void axpy(const ScalarType& alpha, const ThisType& xx)
  {
    assert(has_equal_shape(xx));
    ensure_uniqueness();
    const auto& xx_backend = xx.backend();
    for (size_t ii = 0; ii < backend_->size(); ++ii)
      backend_->operator[](ii) += alpha * xx_backend[ii];
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
    std::fill(yy.begin(), yy.end(), ScalarType(0));
    for (size_t rr = 0; rr < num_rows_; ++rr)
      for (size_t kk = row_pointers_->operator[](rr); kk < row_pointers_->operator[](rr + 1); ++kk)
        yy[rr] += backend_->operator[](kk) * xx[column_indices_->operator[](kk)];
  }

  inline void add_to_entry(const size_t rr, const size_t cc, const ScalarType& value)
  {
    ensure_uniqueness();
    const size_t index = get_index_in_backend(rr, cc);
    assert(index != size_t(-1) && "Entry has to be in the sparsity pattern!");
    backend_->operator[](index) += value;
  }

  inline ScalarType get_entry(const size_t rr, const size_t cc) const
  {
    const size_t index = get_index_in_backend(rr, cc);
    return index == size_t(-1) ? ScalarType(0) : backend_->operator[](index);
  }

  inline void set_entry(const size_t rr, const size_t cc, const ScalarType value)
  {
    ensure_uniqueness();
    const size_t index = get_index_in_backend(rr, cc);
    assert(index != size_t(-1) && "Entry has to be in the sparsity pattern!");
    backend_->operator[](index) = value;
  }

  inline void clear_row(const size_t rr)
  {
    ensure_uniqueness();
    std::fill(backend_->begin() + row_pointers_->operator[](rr),
              backend_->begin() + row_pointers_->operator[](rr + 1),
              ScalarType(0));
  }

  inline void clear_col(const size_t cc)
  {
    ensure_uniqueness();
    for (size_t kk = 0; kk < backend_->size(); ++kk) {
      if (column_indices_->operator[](kk) == cc)
        backend_->operator[](kk) = ScalarType(0);
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
    // iterate over non-zero entries
    for (const auto& entry : *backend_)
      if (DSC::isnan(std::real(entry)) || DSC::isnan(std::imag(entry)) || DSC::isinf(std::abs(entry)))
        return false;
    return true;
  }

  virtual size_t non_zeros() const override final
  {
    return backend_->size();
  }

  virtual SparsityPatternDefault pattern(const bool prune = false,
                                         const typename Common::FloatCmp::DefaultEpsilon<ScalarType>::Type eps =
                                             Common::FloatCmp::DefaultEpsilon<ScalarType>::value()) const override
  {
    SparsityPatternDefault ret(num_rows_);
    for (size_t rr = 0; rr < num_rows_; ++rr) {
      for (size_t kk = row_pointers_->operator[](rr + 1); kk < row_pointers_->operator[](rr + 1); ++kk) {
        const auto& val = backend_->operator[](kk);
        if (!prune || Common::FloatCmp::ne(val, ScalarType(0), eps))
          ret.insert(rr, column_indices_->operator[](kk));
      }
    }
    return ret;
  } // ... pattern(...)

  /// \}

private:
  size_t get_index_in_backend(const size_t rr, const size_t cc) const
  {
    const auto& row_offset = row_pointers_->operator[](rr);
    auto column_indices_iterator = column_indices_->begin();
    column_indices_iterator += row_offset;
    for (size_t kk = row_offset; kk < row_pointers_->operator[](rr + 1); ++kk, ++column_indices_iterator)
      if (*column_indices_iterator == cc)
        return kk;
    return size_t(-1);
  }

  inline void ensure_uniqueness() const
  {
    if (!backend_.unique())
      backend_ = std::make_shared<BackendType>(*backend_);
  } // ... ensure_uniqueness(...)

  size_t num_rows_, num_cols_;
  mutable std::shared_ptr<BackendType> backend_;
  std::shared_ptr<IndexVectorType> row_pointers_;
  std::shared_ptr<IndexVectorType> column_indices_;
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

#endif // DUNE_XT_LA_CONTAINER_COMMON_HH
