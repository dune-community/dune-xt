// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Barbara Verfürth (2015)
//   Felix Schindler  (2014 - 2017)
//   Rene Milk        (2014 - 2016, 2018)
//   Tobias Leibner   (2014 - 2015, 2017 - 2018)

#ifndef DUNE_XT_LA_CONTAINER_EIGEN_SPARSE_HH
#define DUNE_XT_LA_CONTAINER_EIGEN_SPARSE_HH

#include <memory>
#include <type_traits>
#include <vector>
#include <complex>
#include <mutex>

#include <boost/numeric/conversion/cast.hpp>

#if HAVE_EIGEN
#include <dune/xt/common/disable_warnings.hh>
#include <Eigen/SparseCore>
#include <dune/xt/common/reenable_warnings.hh>
#endif

#include <dune/common/typetraits.hh>
#include <dune/common/ftraits.hh>
#include <dune/common/unused.hh>

#include <dune/xt/common/math.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/crtp.hh>
#include <dune/xt/common/float_cmp.hh>

#include "dune/xt/la/container/interfaces.hh"
#include "dune/xt/la/container/pattern.hh"

#include "dense.hh"

namespace Dune {
namespace XT {
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
  static const constexpr Backends backend_type = Backends::eigen_sparse;
  static const constexpr Backends vector_type = Backends::eigen_dense;
  static const constexpr bool sparse = true;
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
  EigenRowMajorSparseMatrix(const size_t rr,
                            const size_t cc,
                            const SparsityPatternDefault& pattern_in,
                            const size_t num_mutexes = 1)
    : backend_(std::make_shared<BackendType>(internal::boost_numeric_cast<EIGEN_size_t>(rr),
                                             internal::boost_numeric_cast<EIGEN_size_t>(cc)))
    , mutexes_(num_mutexes > 0 ? std::make_shared<std::vector<std::mutex>>(num_mutexes) : nullptr)
    , unshareable_(false)
  {
    if (rr > 0 && cc > 0) {
      if (size_t(pattern_in.size()) != rr)
        DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                   "The size of the pattern (" << pattern_in.size() << ") does not match the number of rows of this ("
                                               << rr
                                               << ")!");
      for (size_t row = 0; row < size_t(pattern_in.size()); ++row) {
        backend_->startVec(internal::boost_numeric_cast<EIGEN_size_t>(row));
        const auto& columns = pattern_in.inner(row);
        for (auto& column : columns) {
#ifndef NDEBUG
          if (column >= cc)
            DUNE_THROW(Common::Exceptions::shapes_do_not_match,
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

  explicit EigenRowMajorSparseMatrix(const size_t rr = 0, const size_t cc = 0, const size_t num_mutexes = 1)
    : backend_(std::make_shared<BackendType>(rr, cc))
    , mutexes_(num_mutexes > 0 ? std::make_shared<std::vector<std::mutex>>(num_mutexes) : nullptr)
    , unshareable_(false)
  {
  }

  EigenRowMajorSparseMatrix(const ThisType& other)
    : backend_(other.unshareable_ ? std::make_shared<BackendType>(*other.backend_) : other.backend_)
    , mutexes_(other.unshareable_ ? std::make_shared<std::vector<std::mutex>>(other.mutexes_->size()) : other.mutexes_)
    , unshareable_(false)
  {
  }

  explicit EigenRowMajorSparseMatrix(const BackendType& mat,
                                     const bool prune = false,
                                     const typename Common::FloatCmp::DefaultEpsilon<ScalarType>::Type eps =
                                         Common::FloatCmp::DefaultEpsilon<ScalarType>::value(),
                                     const size_t num_mutexes = 1)
    : mutexes_(num_mutexes > 0 ? std::make_shared<std::vector<std::mutex>>(num_mutexes) : nullptr)
    , unshareable_(false)
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
          if (Common::FloatCmp::ne<Common::FloatCmp::Style::absolute>(val, zero, eps))
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
  explicit EigenRowMajorSparseMatrix(BackendType* backend_ptr, const size_t num_mutexes = 1)
    : backend_(backend_ptr)
    , mutexes_(num_mutexes > 0 ? std::make_shared<std::vector<std::mutex>>(num_mutexes) : nullptr)
    , unshareable_(false)
  {
  }

  explicit EigenRowMajorSparseMatrix(std::shared_ptr<BackendType> backend_ptr, const size_t num_mutexes = 1)
    : backend_(backend_ptr)
    , mutexes_(num_mutexes > 0 ? std::make_shared<std::vector<std::mutex>>(num_mutexes) : nullptr)
    , unshareable_(false)
  {
  }

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      backend_ = other.unshareable_ ? std::make_shared<BackendType>(*other.backend_) : other.backend_;
      mutexes_ =
          other.unshareable_ ? std::make_shared<std::vector<std::mutex>>(other.mutexes_->size()) : other.mutexes_;
      unshareable_ = false;
    }
    return *this;
  }

  /**
   *  \note Does a deep copy.
   */
  ThisType& operator=(const BackendType& other)
  {
    backend_ = std::make_shared<BackendType>(other);
    unshareable_ = false;
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
    return ThisType(*backend_);
  }

  void scal(const ScalarType& alpha)
  {
    auto& backend_ref = backend();
    const internal::VectorLockGuard DUNE_UNUSED(guard)(mutexes_);
    backend_ref *= alpha;
  }

  void axpy(const ScalarType& alpha, const ThisType& xx)
  {
    auto& backend_ref = backend();
    const internal::VectorLockGuard DUNE_UNUSED(guard)(mutexes_);
    if (!has_equal_shape(xx))
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The shape of xx (" << xx.rows() << "x" << xx.cols() << ") does not match the shape of this ("
                                     << rows()
                                     << "x"
                                     << cols()
                                     << ")!");
    backend_ref += alpha * xx.backend();
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
    yy.backend().transpose() = backend() * xx.backend();
  }

  void add_to_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    auto& backend_ref = backend();
    internal::LockGuard DUNE_UNUSED(lock)(mutexes_, ii);
    assert(these_are_valid_indices(ii, jj));
    backend_ref.coeffRef(internal::boost_numeric_cast<EIGEN_size_t>(ii),
                         internal::boost_numeric_cast<EIGEN_size_t>(jj)) += value;
  }

  void set_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    auto& backend_ref = backend();
    assert(these_are_valid_indices(ii, jj));
    backend_ref.coeffRef(internal::boost_numeric_cast<EIGEN_size_t>(ii),
                         internal::boost_numeric_cast<EIGEN_size_t>(jj)) = value;
  }

  ScalarType get_entry(const size_t ii, const size_t jj) const
  {
    assert(ii < rows());
    assert(jj < cols());
    return backend().coeff(internal::boost_numeric_cast<EIGEN_size_t>(ii),
                           internal::boost_numeric_cast<EIGEN_size_t>(jj));
  }

  void clear_row(const size_t ii)
  {
    auto& backend_ref = backend();
    if (ii >= rows())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given ii (" << ii << ") is larger than the rows of this (" << rows() << ")!");
    backend_ref.row(internal::boost_numeric_cast<EIGEN_size_t>(ii)) *= ScalarType(0);
  }

  void clear_col(const size_t jj)
  {
    auto& backend_ref = backend();
    if (jj >= cols())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given jj (" << jj << ") is larger than the cols of this (" << cols() << ")!");
    for (size_t row = 0; internal::boost_numeric_cast<EIGEN_size_t>(row) < backend_ref.outerSize(); ++row) {
      for (typename BackendType::InnerIterator row_it(backend_ref, internal::boost_numeric_cast<EIGEN_size_t>(row));
           row_it;
           ++row_it) {
        const size_t col = row_it.col();
        if (col == jj) {
          backend_ref.coeffRef(internal::boost_numeric_cast<EIGEN_size_t>(row),
                               internal::boost_numeric_cast<EIGEN_size_t>(jj)) = ScalarType(0);
          break;
        } else if (col > jj)
          break;
      }
    }
  } // ... clear_col(...)

  void unit_row(const size_t ii)
  {
    auto& backend_ref = backend();
    if (ii >= cols())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given ii (" << ii << ") is larger than the cols of this (" << cols() << ")!");
    if (ii >= rows())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given ii (" << ii << ") is larger than the rows of this (" << rows() << ")!");
    if (!these_are_valid_indices(ii, ii))
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Diagonal entry (" << ii << ", " << ii << ") is not contained in the sparsity pattern!");
    backend_ref.row(internal::boost_numeric_cast<EIGEN_size_t>(ii)) *= ScalarType(0);
    backend_ref.coeffRef(internal::boost_numeric_cast<EIGEN_size_t>(ii),
                         internal::boost_numeric_cast<EIGEN_size_t>(ii)) = ScalarType(1);
  } // ... unit_row(...)

  void unit_col(const size_t jj)
  {
    auto& backend_ref = backend();
    if (jj >= cols())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given jj (" << jj << ") is larger than the cols of this (" << cols() << ")!");
    if (jj >= rows())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given jj (" << jj << ") is larger than the rows of this (" << rows() << ")!");
    for (size_t row = 0; internal::boost_numeric_cast<EIGEN_size_t>(row) < backend_ref.outerSize(); ++row) {
      for (typename BackendType::InnerIterator row_it(backend_ref, internal::boost_numeric_cast<EIGEN_size_t>(row));
           row_it;
           ++row_it) {
        const size_t col = row_it.col();
        if (col == jj) {
          if (col == row)
            backend_ref.coeffRef(internal::boost_numeric_cast<EIGEN_size_t>(row),
                                 internal::boost_numeric_cast<EIGEN_size_t>(col)) = ScalarType(1);
          else
            backend_ref.coeffRef(internal::boost_numeric_cast<EIGEN_size_t>(row),
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
    for (EIGEN_size_t ii = 0; ii < backend().outerSize(); ++ii) {
      for (InnerIterator it(backend(), ii); it; ++it) {
        if (Common::isnan(it.value()) || Common::isinf(it.value()))
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
      for (EIGEN_size_t row = 0; row < backend().outerSize(); ++row) {
        for (typename BackendType::InnerIterator row_it(backend(), row); row_it; ++row_it) {
          const EIGEN_size_t col = row_it.col();
          const auto val = backend().coeff(row, col);
          if (Common::FloatCmp::ne(val, zero, eps))
            ret.insert(boost::numeric_cast<size_t>(row), boost::numeric_cast<size_t>(col));
        }
      }
    } else {
      for (EIGEN_size_t row = 0; row < backend().outerSize(); ++row) {
        for (typename BackendType::InnerIterator row_it(backend(), row); row_it; ++row_it)
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

  using MatrixInterfaceType::operator+;
  using MatrixInterfaceType::operator-;
  using MatrixInterfaceType::operator+=;
  using MatrixInterfaceType::operator-=;

private:
  bool these_are_valid_indices(const size_t ii, const size_t jj) const
  {
    if (ii >= rows())
      return false;
    if (jj >= cols())
      return false;
    for (size_t row = ii; internal::boost_numeric_cast<EIGEN_size_t>(row) < backend().outerSize(); ++row) {
      for (typename BackendType::InnerIterator row_it(backend(), internal::boost_numeric_cast<EIGEN_size_t>(row));
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

protected:
  inline void ensure_uniqueness() const
  {
    if (!backend_.unique()) {
      assert(!unshareable_);
      const internal::VectorLockGuard DUNE_UNUSED(guard)(mutexes_);
      if (!backend_.unique()) {
        backend_ = std::make_shared<BackendType>(*backend_);
        mutexes_ = std::make_shared<std::vector<std::mutex>>(mutexes_->size());
      }
    }
  } // ... ensure_uniqueness(...)

private:
  mutable std::shared_ptr<BackendType> backend_;
  mutable std::shared_ptr<std::vector<std::mutex>> mutexes_;
  mutable bool unshareable_;
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
} // namespace XT
} // namespace Dune


// begin: this is what we need for the lib
#if DUNE_XT_WITH_PYTHON_BINDINGS && HAVE_EIGEN


extern template class Dune::XT::LA::EigenRowMajorSparseMatrix<double>;
// extern template void Dune::XT::LA::EigenRowMajorSparseMatrix<double>::mv(const
// Dune::XT::LA::EigenDenseVector<double>&,
//                                                                         Dune::XT::LA::EigenDenseVector<double>&)
//                                                                         const;
// extern template void
// Dune::XT::LA::EigenRowMajorSparseMatrix<double>::mv(const Dune::XT::LA::EigenMappedDenseVector<double>&,
//                                                    Dune::XT::LA::EigenMappedDenseVector<double>&) const;
// extern template void
// Dune::XT::LA::EigenRowMajorSparseMatrix<double>::mv(const Dune::XT::LA::EigenDenseVector<double>&,
//                                                    Dune::XT::LA::EigenMappedDenseVector<double>&) const;
// extern template void
// Dune::XT::LA::EigenRowMajorSparseMatrix<double>::mv(const Dune::XT::LA::EigenMappedDenseVector<double>&,
//                                                    Dune::XT::LA::EigenDenseVector<double>&) const;


#endif // DUNE_XT_WITH_PYTHON_BINDINGS && HAVE_EIGEN
// end: this is what we need for the lib


#endif // DUNE_XT_LA_CONTAINER_EIGEN_SPARSE_HH
