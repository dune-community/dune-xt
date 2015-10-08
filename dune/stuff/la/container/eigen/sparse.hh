// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_LA_CONTAINER_EIGEN_SPARSE_HH
#define DUNE_STUFF_LA_CONTAINER_EIGEN_SPARSE_HH

#include <memory>
#include <type_traits>
#include <vector>
#include <complex>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/stuff/common/disable_warnings.hh>
#if HAVE_EIGEN
#include <Eigen/SparseCore>
#endif
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/common/typetraits.hh>
#include <dune/common/ftraits.hh>

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
  typedef typename Dune::FieldTraits<ScalarImp>::field_type ScalarType;
  typedef typename Dune::FieldTraits<ScalarImp>::real_type RealType;
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
  typedef typename Traits::RealType RealType;

private:
  typedef typename BackendType::Index EIGEN_size_t;

public:
  /**
   * \brief This is the constructor of interest which creates a sparse matrix.
   */
  EigenRowMajorSparseMatrix(const size_t rr, const size_t cc, const SparsityPatternDefault& pattern_in)
  {
    backend_ = std::make_shared<BackendType>(internal::boost_numeric_cast<EIGEN_size_t>(rr),
                                             internal::boost_numeric_cast<EIGEN_size_t>(cc));
    if (rr > 0 && cc > 0) {
      if (size_t(pattern_in.size()) != rr)
        DUNE_THROW(Exceptions::shapes_do_not_match,
                   "The size of the pattern (" << pattern_in.size() << ") does not match the number of rows of this ("
                                               << rr
                                               << ")!");
      for (size_t row = 0; row < size_t(pattern_in.size()); ++row) {
        backend_->startVec(internal::boost_numeric_cast<EIGEN_size_t>(row));
        const auto& columns = pattern_in.inner(row);
        for (auto& column : columns) {
#ifndef NDEBUG
          if (column >= cc)
            DUNE_THROW(Exceptions::shapes_do_not_match,
                       "The size of row " << row << " of the pattern does not match the number of columns of this ("
                                          << cc
                                          << ")!");
#endif // NDEBUG
          backend_->insertBackByOuterInner(internal::boost_numeric_cast<EIGEN_size_t>(row),
                                           internal::boost_numeric_cast<EIGEN_size_t>(column));
        }
        // create entry (insertBackByOuterInner() can not handle empty rows)
        if (columns.size() == 0)
          backend_->insertBackByOuterInner(internal::boost_numeric_cast<EIGEN_size_t>(row), 0);
      }
      backend_->finalize();
      backend_->makeCompressed();
    }
  } // EigenRowMajorSparseMatrix(...)

  explicit EigenRowMajorSparseMatrix(const size_t rr = 0, const size_t cc = 0)
  {
    backend_ = std::make_shared<BackendType>(rr, cc);
  }

  /// This constructor is needed for the python bindings.
  explicit EigenRowMajorSparseMatrix(const DUNE_STUFF_SSIZE_T rr, const DUNE_STUFF_SSIZE_T cc = 0)
    : backend_(new BackendType(internal::boost_numeric_cast<EIGEN_size_t>(rr),
                               internal::boost_numeric_cast<EIGEN_size_t>(cc)))
  {
  }

  explicit EigenRowMajorSparseMatrix(const int rr, const int cc = 0)
    : backend_(new BackendType(internal::boost_numeric_cast<EIGEN_size_t>(rr),
                               internal::boost_numeric_cast<EIGEN_size_t>(cc)))
  {
  }

  EigenRowMajorSparseMatrix(const ThisType& other)
    : backend_(other.backend_)
  {
  }

  explicit EigenRowMajorSparseMatrix(const BackendType& mat, const bool prune = false,
                                     const typename Common::FloatCmp::DefaultEpsilon<ScalarType>::Type eps =
                                         Common::FloatCmp::DefaultEpsilon<ScalarType>::value())
  {
    if (prune) {
      // we do this here instead of using pattern(true), since we can build the triplets along the way which is more
      // efficient
      typedef ::Eigen::Triplet<ScalarType> TripletType;
      const ScalarType zero(0);
      std::vector<TripletType> triplets;
      triplets.reserve(mat.nonZeros());
      for (EIGEN_size_t row = 0; row < mat.outerSize(); ++row) {
        for (typename BackendType::InnerIterator row_it(mat, row); row_it; ++row_it) {
          const EIGEN_size_t col = row_it.col();
          const auto val = mat.coeff(row, col);
          if (Stuff::Common::FloatCmp::ne<Stuff::Common::FloatCmp::Style::absolute>(val, zero, eps))
            triplets.emplace_back(row, col, val);
        }
      }
      backend_ = std::make_shared<BackendType>(mat.rows(), mat.cols());
      backend_->setFromTriplets(triplets.begin(), triplets.end());
    } else
      backend_ = std::make_shared<BackendType>(mat);
  } // EigenRowMajorSparseMatrix(...)

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

  template <class T1, class T2>
  inline void mv(const EigenBaseVector<T1, ScalarType>& xx, EigenBaseVector<T2, ScalarType>& yy) const
  {
    yy.backend().transpose() = backend_->operator*(*xx.backend_);
  }

  void add_to_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    assert(these_are_valid_indices(ii, jj));
    backend().coeffRef(internal::boost_numeric_cast<EIGEN_size_t>(ii),
                       internal::boost_numeric_cast<EIGEN_size_t>(jj)) += value;
  }

  void set_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    assert(these_are_valid_indices(ii, jj));
    backend().coeffRef(internal::boost_numeric_cast<EIGEN_size_t>(ii), internal::boost_numeric_cast<EIGEN_size_t>(jj)) =
        value;
  }

  ScalarType get_entry(const size_t ii, const size_t jj) const
  {
    assert(ii < rows());
    assert(jj < cols());
    return backend_->coeff(internal::boost_numeric_cast<EIGEN_size_t>(ii),
                           internal::boost_numeric_cast<EIGEN_size_t>(jj));
  }

  void clear_row(const size_t ii)
  {
    if (ii >= rows())
      DUNE_THROW(Exceptions::index_out_of_range,
                 "Given ii (" << ii << ") is larger than the rows of this (" << rows() << ")!");
    backend().row(internal::boost_numeric_cast<EIGEN_size_t>(ii)) *= ScalarType(0);
  }

  void clear_col(const size_t jj)
  {
    if (jj >= cols())
      DUNE_THROW(Exceptions::index_out_of_range,
                 "Given jj (" << jj << ") is larger than the cols of this (" << cols() << ")!");
    ensure_uniqueness();
    for (size_t row = 0; internal::boost_numeric_cast<EIGEN_size_t>(row) < backend_->outerSize(); ++row) {
      for (typename BackendType::InnerIterator row_it(*backend_, internal::boost_numeric_cast<EIGEN_size_t>(row));
           row_it;
           ++row_it) {
        const size_t col = row_it.col();
        if (col == jj) {
          backend_->coeffRef(internal::boost_numeric_cast<EIGEN_size_t>(row),
                             internal::boost_numeric_cast<EIGEN_size_t>(jj)) = ScalarType(0);
          break;
        } else if (col > jj)
          break;
      }
    }
  } // ... clear_col(...)

  void unit_row(const size_t ii)
  {
    if (ii >= cols())
      DUNE_THROW(Exceptions::index_out_of_range,
                 "Given ii (" << ii << ") is larger than the cols of this (" << cols() << ")!");
    if (ii >= rows())
      DUNE_THROW(Exceptions::index_out_of_range,
                 "Given ii (" << ii << ") is larger than the rows of this (" << rows() << ")!");
    if (!these_are_valid_indices(ii, ii))
      DUNE_THROW(Exceptions::index_out_of_range,
                 "Diagonal entry (" << ii << ", " << ii << ") is not contained in the sparsity pattern!");
    backend().row(internal::boost_numeric_cast<EIGEN_size_t>(ii)) *= ScalarType(0);
    set_entry(ii, ii, ScalarType(1));
  } // ... unit_row(...)

  void unit_col(const size_t jj)
  {
    if (jj >= cols())
      DUNE_THROW(Exceptions::index_out_of_range,
                 "Given jj (" << jj << ") is larger than the cols of this (" << cols() << ")!");
    if (jj >= rows())
      DUNE_THROW(Exceptions::index_out_of_range,
                 "Given jj (" << jj << ") is larger than the rows of this (" << rows() << ")!");
    ensure_uniqueness();
    for (size_t row = 0; internal::boost_numeric_cast<EIGEN_size_t>(row) < backend_->outerSize(); ++row) {
      for (typename BackendType::InnerIterator row_it(*backend_, internal::boost_numeric_cast<EIGEN_size_t>(row));
           row_it;
           ++row_it) {
        const size_t col = row_it.col();
        if (col == jj) {
          if (col == row)
            backend_->coeffRef(internal::boost_numeric_cast<EIGEN_size_t>(row),
                               internal::boost_numeric_cast<EIGEN_size_t>(col)) = ScalarType(1);
          else
            backend_->coeffRef(internal::boost_numeric_cast<EIGEN_size_t>(row),
                               internal::boost_numeric_cast<EIGEN_size_t>(jj)) = ScalarType(0);
          break;
        } else if (col > jj)
          break;
      }
    }
  } // ... unit_col(...)

  bool valid() const
  {
    // iterate over non-zero entries
    typedef typename BackendType::InnerIterator InnerIterator;
    for (EIGEN_size_t ii = 0; ii < backend_->outerSize(); ++ii) {
      for (InnerIterator it(*backend_, ii); it; ++it) {
        if (DSC::isnan(std::real(it.value())) || DSC::isnan(std::imag(it.value())) || DSC::isinf(std::abs(it.value())))
          return false;
      }
    }
    return true;
  }

  virtual size_t non_zeros() const override final
  {
    return backend_->nonZeros();
  }

  virtual SparsityPatternDefault
  pattern(const bool prune = false,
          const ScalarType eps = Common::FloatCmp::DefaultEpsilon<ScalarType>::value()) const override
  {
    SparsityPatternDefault ret(rows());
    const auto zero = typename Common::FloatCmp::DefaultEpsilon<ScalarType>::Type(0);
    if (prune) {
      for (EIGEN_size_t row = 0; row < backend_->outerSize(); ++row) {
        for (typename BackendType::InnerIterator row_it(*backend_, row); row_it; ++row_it) {
          const EIGEN_size_t col = row_it.col();
          const auto val = backend_->coeff(row, col);
          if (Common::FloatCmp::ne(val, zero, eps))
            ret.insert(boost::numeric_cast<size_t>(row), boost::numeric_cast<size_t>(col));
        }
      }
    } else {
      for (EIGEN_size_t row = 0; row < backend_->outerSize(); ++row) {
        for (typename BackendType::InnerIterator row_it(*backend_, row); row_it; ++row_it)
          ret.insert(boost::numeric_cast<size_t>(row), boost::numeric_cast<size_t>(row_it.col()));
      }
    }
    ret.sort();
    return ret;
  } // ... pattern(...)

  virtual ThisType
  pruned(const ScalarType eps = Common::FloatCmp::DefaultEpsilon<ScalarType>::value()) const override final
  {
    return ThisType(*backend_, true, eps);
  }

  /// \}

private:
  bool these_are_valid_indices(const size_t ii, const size_t jj) const
  {
    if (ii >= rows())
      return false;
    if (jj >= cols())
      return false;
    for (size_t row = ii; internal::boost_numeric_cast<EIGEN_size_t>(row) < backend_->outerSize(); ++row) {
      for (typename BackendType::InnerIterator row_it(*backend_, internal::boost_numeric_cast<EIGEN_size_t>(row));
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
  static_assert(AlwaysFalse<ScalarImp>::value, "You are missing Eigen!");
};

#endif // HAVE_EIGEN

} // namespace LA
namespace Common {

#if HAVE_EIGEN

template <class T>
struct MatrixAbstraction<LA::EigenRowMajorSparseMatrix<T>>
    : public LA::internal::MatrixAbstractionBase<LA::EigenRowMajorSparseMatrix<T>>
{
};

#endif // HAVE_EIGEN

} // namespace Common
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_LA_CONTAINER_EIGEN_SPARSE_HH
