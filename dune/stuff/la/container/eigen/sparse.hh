// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_LA_CONTAINER_EIGEN_SPARSE_HH
#define DUNE_STUFF_LA_CONTAINER_EIGEN_SPARSE_HH

#include <memory>
#include <type_traits>
#include <vector>

#if HAVE_EIGEN
#include <dune/stuff/common/disable_warnings.hh>
#include <Eigen/SparseCore>
#include <dune/stuff/common/reenable_warnings.hh>
#endif // HAVE_EIGEN

#include <dune/common/typetraits.hh>

#include <dune/stuff/aliases.hh>
#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/crtp.hh>
#include <dune/stuff/common/float_cmp.hh>

#include "dune/stuff/la/container/interfaces.hh"
#include "dune/stuff/la/container/pattern.hh"

#include "dense.hh"

namespace Dune {
namespace Stuff {
namespace LA {


// forwards
template <class ScalarType>
class EigenRowMajorSparseMatrix;

class EigenMatrixInterfaceDynamic
{
};


#if HAVE_EIGEN


namespace internal {


/**
 * \brief Traits for EigenRowMajorSparseMatrix.
 */
template <class ScalarImp = double>
class EigenRowMajorSparseMatrixTraits
{
public:
  typedef ScalarImp ScalarType;
  typedef EigenRowMajorSparseMatrix<ScalarType> derived_type;
  typedef typename ::Eigen::SparseMatrix<ScalarType, ::Eigen::RowMajor> BackendType;
}; // class RowMajorSparseMatrixTraits


} // namespace internal


/**
 * \brief A sparse matrix implementation of the MatrixInterface with row major memory layout.
 */
template <class ScalarImp = double>
class EigenRowMajorSparseMatrix
    : public MatrixInterface<internal::EigenRowMajorSparseMatrixTraits<ScalarImp>, ScalarImp>,
      public ProvidesBackend<internal::EigenRowMajorSparseMatrixTraits<ScalarImp>>
{
  typedef EigenRowMajorSparseMatrix<ScalarImp> ThisType;
  typedef MatrixInterface<internal::EigenRowMajorSparseMatrixTraits<ScalarImp>, ScalarImp> MatrixInterfaceType;
  static_assert(!std::is_same<DUNE_STUFF_SSIZE_T, int>::value,
                "You have to manually disable the constructor below which uses DUNE_STUFF_SSIZE_T!");

public:
  typedef internal::EigenRowMajorSparseMatrixTraits<ScalarImp> Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::ScalarType ScalarType;

  /**
   * \brief This is the constructor of interest which creates a sparse matrix.
   */
  EigenRowMajorSparseMatrix(const size_t rr, const size_t cc, const SparsityPatternDefault& pattern)
  {
    backend_ = std::make_shared<BackendType>(assert_is_IndexType_compatible_and_convert(rr),
                                             assert_is_IndexType_compatible_and_convert(cc));
    if (rr > 0 && cc > 0) {
      if (size_t(pattern.size()) != rr)
        DUNE_THROW(Exceptions::shapes_do_not_match,
                   "The size of the pattern (" << pattern.size() << ") does not match the number of rows of this ("
                                               << rr
                                               << ")!");
      for (size_t row = 0; row < size_t(pattern.size()); ++row) {
        backend_->startVec(assert_is_IndexType_compatible_and_convert(row));
        const auto& columns = pattern.inner(row);
        for (auto& column : columns) {
#ifndef NDEBUG
          if (column >= cc)
            DUNE_THROW(Exceptions::shapes_do_not_match,
                       "The size of row " << row << " of the pattern does not match the number of columns of this ("
                                          << cc
                                          << ")!");
#endif // NDEBUG
          backend_->insertBackByOuterInner(assert_is_IndexType_compatible_and_convert(row),
                                           assert_is_IndexType_compatible_and_convert(column));
        }
        // create entry (insertBackByOuterInner() can not handle empty rows)
        if (columns.size() == 0)
          backend_->insertBackByOuterInner(assert_is_IndexType_compatible_and_convert(row), 0);
      }
      backend_->finalize();
      backend_->makeCompressed();
    }
  }

  explicit EigenRowMajorSparseMatrix(const size_t rr = 0, const size_t cc = 0)
  {
    backend_ = std::make_shared<BackendType>(assert_is_IndexType_compatible_and_convert(rr),
                                             assert_is_IndexType_compatible_and_convert(cc));
  }

  /// This constructor is needed for the python bindings.
  explicit EigenRowMajorSparseMatrix(const DUNE_STUFF_SSIZE_T rr, const DUNE_STUFF_SSIZE_T cc = 0)
    : backend_(new BackendType(MatrixInterfaceType::assert_is_size_t_compatible_and_convert(rr),
                               MatrixInterfaceType::assert_is_size_t_compatible_and_convert(cc)))
  {
  }

  explicit EigenRowMajorSparseMatrix(const int rr, const int cc = 0)
    : EigenRowMajorSparseMatrix(MatrixInterfaceType::assert_is_size_t_compatible_and_convert(rr),
                                MatrixInterfaceType::assert_is_size_t_compatible_and_convert(cc))
  {
    //    backend_ = std::make_shared<BackendType>(assert_is_IndexType_compatible_and_convert(rr),
    //                                             assert_is_IndexType_compatible_and_convert(cc));
  }

  EigenRowMajorSparseMatrix(const ThisType& other)
    : backend_(other.backend_)
  {
  }

  explicit EigenRowMajorSparseMatrix(const BackendType& other)
    : backend_(new BackendType(other))
  {
  }

  /**
   *  \note Takes ownership of backend_ptr in the sense that you must not delete it afterwards!
   */
  explicit EigenRowMajorSparseMatrix(BackendType* backend_ptr)
    : backend_(backend_ptr)
  {
  }

  explicit EigenRowMajorSparseMatrix(std::shared_ptr<BackendType> backend_ptr)
    : backend_(backend_ptr)
  {
  }

  ThisType& operator=(const ThisType& other)
  {
    backend_ = other.backend_;
    return *this;
  } // ... operator=(...)

  /**
   *  \note Does a deep copy.
   */
  ThisType& operator=(const BackendType& other)
  {
    backend_ = std::make_shared<BackendType>(other);
    return *this;
  } // ... operator=(...)

  /**
   * \defgroup backend ´´These methods are required by the ProvidesBackend interface.``
   * \{
   */

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
  /**
   * \}
   */

  /**
   * \defgroup container ´´These methods are required by ContainerInterface.``
   * \{
   */

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
    if (!has_equal_shape(xx))
      DUNE_THROW(Exceptions::shapes_do_not_match,
                 "The shape of xx (" << xx.rows() << "x" << xx.cols() << ") does not match the shape of this ("
                                     << rows()
                                     << "x"
                                     << cols()
                                     << ")!");
    const auto& xx_ref = *(xx.backend_);
    backend() += alpha * xx_ref;
  } // ... axpy(...)

  bool has_equal_shape(const ThisType& other) const
  {
    return (rows() == other.rows()) && (cols() == other.cols());
  }
  /**
   * \}
   */

  /**
   * \defgroup matrix_required ´´These methods are required by MatrixInterface.``
   * \{
   */

  inline size_t rows() const
  {
    return backend_->rows();
  }

  inline size_t cols() const
  {
    return backend_->cols();
  }

  template <class T1, class T2>
  inline void mv(const EigenBaseVector<T1, ScalarType>& xx, EigenBaseVector<T2, ScalarType>& yy) const
  {
    yy.backend().transpose() = backend_->operator*(*xx.backend_);
  }

  void add_to_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    assert(these_are_valid_indices(ii, jj));
    backend().coeffRef(assert_is_IndexType_compatible_and_convert(ii),
                       assert_is_IndexType_compatible_and_convert(jj)) += value;
  } // ... add_to_entry(...)

  void set_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    assert(these_are_valid_indices(ii, jj));
    backend().coeffRef(assert_is_IndexType_compatible_and_convert(ii), assert_is_IndexType_compatible_and_convert(jj)) =
        value;
  } // ... set_entry(...)

  ScalarType get_entry(const size_t ii, const size_t jj) const
  {
    assert(ii < rows());
    assert(jj < cols());
    return backend_->coeff(assert_is_IndexType_compatible_and_convert(ii),
                           assert_is_IndexType_compatible_and_convert(jj));
  } // ... get_entry(...)

  void clear_row(const size_t ii)
  {
    if (ii >= rows())
      DUNE_THROW(Exceptions::index_out_of_range,
                 "Given ii (" << ii << ") is larger than the rows of this (" << rows() << ")!");
    backend().row(assert_is_IndexType_compatible_and_convert(ii)) *= ScalarType(0);
  } // ... clear_row(...)

  void clear_col(const size_t jj)
  {
    if (jj >= cols())
      DUNE_THROW(Exceptions::index_out_of_range,
                 "Given jj (" << jj << ") is larger than the cols of this (" << cols() << ")!");
    ensure_uniqueness();
    for (size_t row = 0; assert_is_IndexType_compatible_and_convert(row) < backend_->outerSize(); ++row) {
      for (typename BackendType::InnerIterator row_it(*backend_, assert_is_IndexType_compatible_and_convert(row));
           row_it;
           ++row_it) {
        const size_t col = row_it.col();
        if (col == jj) {
          backend_->coeffRef(assert_is_IndexType_compatible_and_convert(row),
                             assert_is_IndexType_compatible_and_convert(jj)) = ScalarType(0);
          break;
        } else if (col > jj)
          break;
      }
    }
  } // ... clear_col(...)

  void unit_row(const size_t ii)
  {
    if (ii >= rows())
      DUNE_THROW(Exceptions::index_out_of_range,
                 "Given ii (" << ii << ") is larger than the rows of this (" << rows() << ")!");
    if (!these_are_valid_indices(ii, ii))
      DUNE_THROW(Exceptions::index_out_of_range,
                 "Diagonal entry (" << ii << ", " << ii << ") is not contained in the sparsity pattern!");
    backend().row(assert_is_IndexType_compatible_and_convert(ii)) *= ScalarType(0);
    set_entry(ii, ii, ScalarType(1));
  } // ... unit_row(...)

  void unit_col(const size_t jj)
  {
    if (jj >= rows())
      DUNE_THROW(Exceptions::index_out_of_range,
                 "Given jj (" << jj << ") is larger than the cols of this (" << cols() << ")!");
    ensure_uniqueness();
    for (size_t row = 0; assert_is_IndexType_compatible_and_convert(row) < backend_->outerSize(); ++row) {
      for (typename BackendType::InnerIterator row_it(*backend_, assert_is_IndexType_compatible_and_convert(row));
           row_it;
           ++row_it) {
        const size_t col = row_it.col();
        if (col == jj) {
          if (col == row)
            backend_->coeffRef(assert_is_IndexType_compatible_and_convert(row),
                               assert_is_IndexType_compatible_and_convert(col)) = ScalarType(1);
          else
            backend_->coeffRef(assert_is_IndexType_compatible_and_convert(row),
                               assert_is_IndexType_compatible_and_convert(jj)) = ScalarType(0);
          break;
        } else if (col > jj)
          break;
      }
    }
  } // ... unit_col(...)

  bool valid() const
  {
    // serialize matrix (no copy done here)
    auto& non_const_ref = const_cast<BackendType&>(*backend_);
    return EigenMappedDenseVector<ScalarType>(non_const_ref.valuePtr(), non_const_ref.nonZeros()).valid();
  } // ... valid(...)

  /**
   * \}
   */

private:
  typedef typename BackendType::Index IndexType;

  bool these_are_valid_indices(const size_t ii, const size_t jj) const
  {
    if (ii >= rows())
      return false;
    if (jj >= cols())
      return false;
    for (size_t row = ii; assert_is_IndexType_compatible_and_convert(row) < backend_->outerSize(); ++row) {
      for (typename BackendType::InnerIterator row_it(*backend_, assert_is_IndexType_compatible_and_convert(row));
           row_it;
           ++row_it) {
        const size_t col = row_it.col();
        if ((ii == row) && (jj == col))
          return true;
        else if ((row > ii) && (col > jj))
          return false;
      }
    }
    return false;
  } // ... these_are_valid_indices(...)

  template <class SizeType>
  static inline IndexType assert_is_IndexType_compatible_and_convert(const SizeType& size)
  {
    assert(size < std::numeric_limits<IndexType>::max());
    return IndexType(size);
  } // ... assert_is_IndexType_compatible_and_convert(...)

  inline void ensure_uniqueness() const
  {
    if (!backend_.unique())
      backend_ = std::make_shared<BackendType>(*backend_);
  } // ... ensure_uniqueness(...)

  mutable std::shared_ptr<BackendType> backend_;
}; // class EigenRowMajorSparseMatrix


#else // HAVE_EIGEN


template <class ScalarImp>
class EigenRowMajorSparseMatrix
{
  static_assert(Dune::AlwaysFalse<ScalarImp>::value, "You are missing Eigen!");
};


#endif // HAVE_EIGEN

} // namespace LA
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_LA_CONTAINER_EIGEN_SPARSE_HH
