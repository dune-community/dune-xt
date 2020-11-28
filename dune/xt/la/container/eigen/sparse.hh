// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Barbara Verfürth (2015)
//   Felix Schindler  (2014 - 2018, 2020)
//   René Fritze      (2014 - 2016, 2018 - 2020)
//   Tobias Leibner   (2014 - 2015, 2017 - 2020)

#ifndef DUNE_XT_LA_CONTAINER_EIGEN_SPARSE_HH
#define DUNE_XT_LA_CONTAINER_EIGEN_SPARSE_HH

#include <memory>
#include <type_traits>
#include <vector>
#include <complex>
#include <mutex>

#include <boost/numeric/conversion/cast.hpp>

#if HAVE_EIGEN
#  include <dune/xt/common/disable_warnings.hh>
#  include <Eigen/SparseCore>
#  include <dune/xt/common/reenable_warnings.hh>
#endif

#include <dune/common/typetraits.hh>
#include <dune/common/ftraits.hh>

#include <dune/xt/common/crtp.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/math.hh>
#include <dune/xt/common/numeric_cast.hh>

#include "dune/xt/la/container/interfaces.hh"
#include "dune/xt/la/container/pattern.hh"

#include "dense.hh"

namespace Dune {
namespace XT {
namespace LA {


// forwards
template <class ScalarType>
class EigenRowMajorSparseMatrix;

template <class ScalarType>
class EigenColMajorSparseMatrix;

class EigenMatrixInterfaceDynamic
{};


#if HAVE_EIGEN

namespace internal {


/**
 * \brief Traits for EigenRowMajorSparseMatrix.
 */
template <class ScalarImp>
struct EigenRowMajorSparseMatrixTraits
  : public MatrixTraitsBase<ScalarImp,
                            EigenRowMajorSparseMatrix<ScalarImp>,
                            ::Eigen::SparseMatrix<ScalarImp, ::Eigen::RowMajor>,
                            Backends::eigen_sparse,
                            Backends::eigen_dense,
                            true>
{
  static constexpr ::Eigen::StorageOptions storage_scheme = ::Eigen::RowMajor;
};

/**
 * \brief Traits for EigenRowMajorSparseMatrix.
 */
template <class ScalarImp>
struct EigenColMajorSparseMatrixTraits
  : public MatrixTraitsBase<ScalarImp,
                            EigenColMajorSparseMatrix<ScalarImp>,
                            ::Eigen::SparseMatrix<ScalarImp, ::Eigen::ColMajor>,
                            Backends::eigen_sparse,
                            Backends::eigen_dense,
                            true>
{
  static constexpr ::Eigen::StorageOptions storage_scheme = ::Eigen::ColMajor;
};

/**
 * \brief A sparse matrix implementation of the MatrixInterface
 */
template <class TraitsImp, class ScalarImp = typename TraitsImp::ScalarType>
class EigenSparseMatrixBase
  : public MatrixInterface<TraitsImp, ScalarImp>
  , public ProvidesBackend<TraitsImp>
{
  using ThisType = EigenSparseMatrixBase;
  using InterfaceType = MatrixInterface<TraitsImp, ScalarImp>;

public:
  using DerivedType = typename InterfaceType::derived_type;
  using Traits = typename InterfaceType::Traits;
  using typename InterfaceType::RealType;
  using typename InterfaceType::ScalarType;
  using typename ProvidesBackend<TraitsImp>::BackendType;
  static constexpr ::Eigen::StorageOptions storage_scheme = TraitsImp::storage_scheme;

private:
  using MutexesType = typename TraitsImp::MutexesType;
  using EIGEN_size_t = typename BackendType::Index;

public:
  /**
   * \brief This is the constructor of interest which creates a sparse matrix.
   */
  EigenSparseMatrixBase(const size_t rr,
                        const size_t cc,
                        const SparsityPatternDefault& pattern_in,
                        const size_t num_mutexes = 1)
    : backend_(
        std::make_shared<BackendType>(Common::numeric_cast<EIGEN_size_t>(rr), Common::numeric_cast<EIGEN_size_t>(cc)))
    , mutexes_(std::make_unique<MutexesType>(num_mutexes))
  {
    if (rr > 0 && cc > 0) {
      if (size_t(pattern_in.size()) != rr)
        DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                   "The size of the pattern (" << pattern_in.size() << ") does not match the number of rows of this ("
                                               << rr << ")!");

      using TripletType = ::Eigen::Triplet<ScalarType>;
      std::vector<TripletType> triplets;
      triplets.reserve(pattern_in.size() * 9);
      for (size_t row = 0; row < rr; ++row) {
        for (const size_t col : pattern_in.inner(row)) {
#  ifndef NDEBUG
          if (col >= cc)
            DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                       "The size of row " << row << " of the pattern does not match the number of columns of this ("
                                          << cc << ")!");
#  endif // NDEBUG
          triplets.emplace_back(static_cast<EIGEN_size_t>(row), static_cast<EIGEN_size_t>(col), ScalarType(0));
        }
      }
      backend_ = std::make_shared<BackendType>(rr, cc);
      backend_->setFromTriplets(triplets.begin(), triplets.end());
      backend_->makeCompressed();
    }
  } // EigenSparseMatrixBase(...)

  explicit EigenSparseMatrixBase(const size_t rr = 0, const size_t cc = 0, const size_t num_mutexes = 1)
    : backend_(
        std::make_shared<BackendType>(Common::numeric_cast<EIGEN_size_t>(rr), Common::numeric_cast<EIGEN_size_t>(cc)))
    , mutexes_(std::make_unique<MutexesType>(num_mutexes))
  {}

  explicit EigenSparseMatrixBase(const size_t rr, const size_t cc, const ScalarType& val, const size_t num_mutexes = 1)
    : backend_(
        std::make_shared<BackendType>(Common::numeric_cast<EIGEN_size_t>(rr), Common::numeric_cast<EIGEN_size_t>(cc)))
    , mutexes_(std::make_unique<MutexesType>(num_mutexes))
  {
    if (val != 0.) {
      for (size_t ii = 0; ii < rows(); ++ii)
        for (size_t jj = 0; jj < cols(); ++jj) {
          backend_->insert(ii, jj);
          set_entry(ii, jj, val);
        }
    }
  }

  EigenSparseMatrixBase(const ThisType& other)
    : backend_(std::make_shared<BackendType>(*other.backend_))
    , mutexes_(std::make_unique<MutexesType>(other.mutexes_->size()))
  {}

  explicit EigenSparseMatrixBase(const BackendType& mat,
                                 const bool prune = false,
                                 const typename Common::FloatCmp::DefaultEpsilon<ScalarType>::Type eps =
                                     Common::FloatCmp::DefaultEpsilon<ScalarType>::value(),
                                 const size_t num_mutexes = 1)
    : mutexes_(std::make_unique<MutexesType>(num_mutexes))
  {
    if (prune) {
      // we do this here instead of using pattern(true), since we can build the triplets along the way which is more
      // efficient
      using TripletType = ::Eigen::Triplet<ScalarType>;
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
  } // EigenSparseMatrixBase(...)

  /**
   *  \note Takes ownership of backend_ptr in the sense that you must not delete it afterwards!
   */
  explicit EigenSparseMatrixBase(BackendType* backend_ptr, const size_t num_mutexes = 1)
    : backend_(backend_ptr)
    , mutexes_(std::make_unique<MutexesType>(num_mutexes))
  {}

  explicit EigenSparseMatrixBase(std::shared_ptr<BackendType> backend_ptr, const size_t num_mutexes = 1)
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
    return *backend_;
  }

  const BackendType& backend() const
  {
    return *backend_;
  }

  /// \}
  /// \name Required by ContainerInterface.
  /// \{

  DerivedType copy() const
  {
    return DerivedType(*backend_);
  }

  void scal(const ScalarType& alpha)
  {
    [[maybe_unused]] const internal::VectorLockGuard guard(*mutexes_);
    backend() *= alpha;
  }

  void set_to_zero()
  {
    backend() *= 0.;
  }

  template <class T>
  void axpy(const ScalarType& alpha, const EigenSparseMatrixBase<T, ScalarImp>& xx)
  {
    if (!has_equal_shape(xx))
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The shape of xx (" << xx.rows() << "x" << xx.cols() << ") does not match the shape of this ("
                                     << rows() << "x" << cols() << ")!");
    [[maybe_unused]] const internal::VectorLockGuard guard(*mutexes_);
    backend() += alpha * xx.backend();
  } // ... axpy(...)

  template <class T>
  bool has_equal_shape(const MatrixInterface<T, ScalarImp>& other) const
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

  template <class V1, class V2>
  inline std::enable_if_t<XT::Common::is_vector<V1>::value && XT::Common::is_vector<V2>::value
                              && (!is_eigen_vector<V1>::value || !is_eigen_vector<V2>::value),
                          void>
  mv(const V1& xx, V2& yy) const
  {
    EigenDenseVector<ScalarType> xx_eigen(xx.size()), yy_eigen(yy.size());
    for (size_t ii = 0; ii < xx.size(); ++ii)
      xx_eigen.set_entry(ii, XT::Common::VectorAbstraction<V1>::get_entry(xx, ii));
    mv(xx_eigen, yy_eigen);
    for (size_t ii = 0; ii < yy.size(); ++ii)
      XT::Common::VectorAbstraction<V2>::set_entry(yy, ii, yy_eigen[ii]);
  }

  template <class T1, class T2>
  inline void mtv(const EigenBaseVector<T1, ScalarType>& xx, EigenBaseVector<T2, ScalarType>& yy) const
  {
    yy.backend().transpose() = backend().transpose() * xx.backend();
  }

  template <class V1, class V2>
  inline std::enable_if_t<XT::Common::is_vector<V1>::value && XT::Common::is_vector<V2>::value
                              && !is_eigen_vector<V1>::value,
                          void>
  mtv(const V1& xx, V2& yy) const
  {
    EigenDenseVector<ScalarType> xx_eigen(xx.size()), yy_eigen(yy.size());
    for (size_t ii = 0; ii < xx.size(); ++ii)
      xx_eigen.set_entry(ii, XT::Common::VectorAbstraction<V1>::get_entry(xx, ii));
    mtv(xx_eigen, yy_eigen);
    for (size_t ii = 0; ii < yy.size(); ++ii)
      XT::Common::VectorAbstraction<V2>::set_entry(yy, ii, yy_eigen[ii]);
  }

  void add_to_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    assert(these_are_valid_indices(ii, jj));
    [[maybe_unused]] const internal::LockGuard lock(*mutexes_, ii, rows());
    backend().coeffRef(static_cast<EIGEN_size_t>(ii), static_cast<EIGEN_size_t>(jj)) += value;
  }

  void set_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    assert(these_are_valid_indices(ii, jj));
    backend().coeffRef(static_cast<EIGEN_size_t>(ii), static_cast<EIGEN_size_t>(jj)) = value;
  }

  ScalarType get_entry(const size_t ii, const size_t jj) const
  {
    assert(ii < rows());
    assert(jj < cols());
    return backend().coeff(static_cast<EIGEN_size_t>(ii), static_cast<EIGEN_size_t>(jj));
  }

  void clear_row(const size_t ii)
  {
    if (ii >= rows())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given ii (" << ii << ") is larger than the rows of this (" << rows() << ")!");
    if constexpr (storage_scheme == ::Eigen::RowMajor) {
      backend().row(static_cast<EIGEN_size_t>(ii)) *= ScalarType(0);
    } else {
      for (size_t col = 0; static_cast<EIGEN_size_t>(col) < backend().outerSize(); ++col) {
        for (typename BackendType::InnerIterator col_it(backend(), static_cast<EIGEN_size_t>(col)); col_it; ++col_it) {
          const size_t row = col_it.row();
          if (row == ii) {
            backend().coeffRef(static_cast<EIGEN_size_t>(row), static_cast<EIGEN_size_t>(col)) = ScalarType(0);
            break;
          } else if (row > ii)
            break;
        }
      }
    }
  }

  void clear_col(const size_t jj)
  {
    if (jj >= cols())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given jj (" << jj << ") is larger than the cols of this (" << cols() << ")!");
    if constexpr (storage_scheme == ::Eigen::RowMajor) {
      for (size_t row = 0; static_cast<EIGEN_size_t>(row) < backend().outerSize(); ++row) {
        for (typename BackendType::InnerIterator row_it(backend(), static_cast<EIGEN_size_t>(row)); row_it; ++row_it) {
          const size_t col = row_it.col();
          if (col == jj) {
            backend().coeffRef(static_cast<EIGEN_size_t>(row), static_cast<EIGEN_size_t>(jj)) = ScalarType(0);
            break;
          } else if (col > jj)
            break;
        }
      }
    } else {
      backend().col(static_cast<EIGEN_size_t>(jj)) *= ScalarType(0);
    }
  } // ... clear_col(...)

  void unit_row(const size_t ii)
  {
    if (!these_are_valid_indices(ii, ii))
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Diagonal entry (" << ii << ", " << ii << ") is not contained in the sparsity pattern!");
    clear_row(ii);
    backend().coeffRef(static_cast<EIGEN_size_t>(ii), static_cast<EIGEN_size_t>(ii)) = ScalarType(1);
  } // ... unit_row(...)

  void unit_col(const size_t jj)
  {
    if (!these_are_valid_indices(jj, jj))
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Diagonal entry (" << jj << ", " << jj << ") is not contained in the sparsity pattern!");
    clear_col(jj);
    backend().coeffRef(static_cast<EIGEN_size_t>(jj), static_cast<EIGEN_size_t>(jj)) = ScalarType(1);
  } // ... unit_col(...)

  bool valid() const
  {
    // iterate over non-zero entries
    using InnerIterator = typename BackendType::InnerIterator;
    for (EIGEN_size_t ii = 0; ii < backend().outerSize(); ++ii) {
      for (InnerIterator it(backend(), ii); it; ++it) {
        if (Common::isnan(it.value()) || Common::isinf(it.value()))
          return false;
      }
    }
    return true;
  }

  size_t non_zeros() const override
  {
    return backend_->nonZeros();
  }

  SparsityPatternDefault
  pattern(const bool prune = false,
          const ScalarType eps = Common::FloatCmp::DefaultEpsilon<ScalarType>::value()) const override
  {
    SparsityPatternDefault ret(rows());
    const auto zero = typename Common::FloatCmp::DefaultEpsilon<ScalarType>::Type(0);
    if (prune) {
      if constexpr (storage_scheme == ::Eigen::RowMajor) {
        for (EIGEN_size_t row = 0; row < backend().outerSize(); ++row) {
          for (typename BackendType::InnerIterator row_it(backend(), row); row_it; ++row_it) {
            const EIGEN_size_t col = row_it.col();
            const auto val = backend().coeff(row, col);
            if (Common::FloatCmp::ne(val, zero, eps))
              ret.insert(static_cast<size_t>(row), static_cast<size_t>(col));
          }
        }
      } else {
        for (EIGEN_size_t col = 0; col < backend().outerSize(); ++col) {
          for (typename BackendType::InnerIterator col_it(backend(), col); col_it; ++col_it) {
            const EIGEN_size_t row = col_it.row();
            const auto val = backend().coeff(row, col);
            if (Common::FloatCmp::ne(val, zero, eps))
              ret.insert(static_cast<size_t>(row), static_cast<size_t>(col));
          }
        }
      }
    } else {
      if constexpr (storage_scheme == ::Eigen::RowMajor) {
        for (EIGEN_size_t row = 0; row < backend().outerSize(); ++row) {
          for (typename BackendType::InnerIterator row_it(backend(), row); row_it; ++row_it)
            ret.insert(static_cast<size_t>(row), static_cast<size_t>(row_it.col()));
        }
      } else {
        for (EIGEN_size_t col = 0; col < backend().outerSize(); ++col) {
          for (typename BackendType::InnerIterator col_it(backend(), col); col_it; ++col_it)
            ret.insert(static_cast<size_t>(col_it.row()), static_cast<size_t>(col));
        }
      }
    }
    ret.sort();
    return ret;
  } // ... pattern(...)

  DerivedType pruned(const ScalarType eps = Common::FloatCmp::DefaultEpsilon<ScalarType>::value()) const override
  {
    return DerivedType(*backend_, true, eps);
  }

  /// \}

  ScalarType* entries()
  {
    backend_->makeCompressed();
    return backend().valuePtr();
  }

  const ScalarType* entries() const
  {
    backend_->makeCompressed();
    return backend().valuePtr();
  }

  int* outer_index_ptr()
  {
    backend_->makeCompressed();
    return backend().outerIndexPtr();
  }

  const int* outer_index_ptr() const
  {
    return backend().outerIndexPtr();
    backend_->makeCompressed();
  }

  int* inner_index_ptr()
  {
    backend_->makeCompressed();
    return backend().innerIndexPtr();
  }

  const int* inner_index_ptr() const
  {
    backend_->makeCompressed();
    return backend().innerIndexPtr();
  }

  void set_num_mutexes(const size_t num_mutexes)
  {
    mutexes_ = std::make_unique<MutexesType>(num_mutexes);
  }

  using InterfaceType::operator+;
  using InterfaceType::operator-;
  using InterfaceType::operator+=;
  using InterfaceType::operator-=;

private:
  bool these_are_valid_indices(const size_t ii, const size_t jj) const
  {
    if (ii >= rows())
      return false;
    if (jj >= cols())
      return false;
    if constexpr (storage_scheme == ::Eigen::RowMajor) {
      for (size_t row = ii; static_cast<EIGEN_size_t>(row) < backend().outerSize(); ++row) {
        for (typename BackendType::InnerIterator row_it(backend(), static_cast<EIGEN_size_t>(row)); row_it; ++row_it) {
          const size_t col = row_it.col();
          if ((ii == row) && (jj == col))
            return true;
          else if ((row > ii) && (col > jj))
            return false;
        }
      }
    } else {
      for (size_t col = jj; static_cast<EIGEN_size_t>(col) < backend().outerSize(); ++col) {
        for (typename BackendType::InnerIterator col_it(backend(), static_cast<EIGEN_size_t>(col)); col_it; ++col_it) {
          const size_t row = col_it.row();
          if ((ii == row) && (jj == col))
            return true;
          else if ((row > ii) && (col > jj))
            return false;
        }
      }
    }
    return false;
  } // ... these_are_valid_indices(...)

  std::shared_ptr<BackendType> backend_;
  std::unique_ptr<MutexesType> mutexes_;
}; // class EigenSparseMatrixBase


} // namespace internal


template <class ScalarImp>
class EigenRowMajorSparseMatrix
  : public internal::EigenSparseMatrixBase<internal::EigenRowMajorSparseMatrixTraits<ScalarImp>, ScalarImp>
{
  using BaseType = internal::EigenSparseMatrixBase<internal::EigenRowMajorSparseMatrixTraits<ScalarImp>, ScalarImp>;

public:
  template <class... Args>
  EigenRowMajorSparseMatrix(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {}
};

template <class ScalarImp>
class EigenColMajorSparseMatrix
  : public internal::EigenSparseMatrixBase<internal::EigenColMajorSparseMatrixTraits<ScalarImp>, ScalarImp>
{
  using BaseType = internal::EigenSparseMatrixBase<internal::EigenColMajorSparseMatrixTraits<ScalarImp>, ScalarImp>;

public:
  template <class... Args>
  EigenColMajorSparseMatrix(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {}
};


#else // HAVE_EIGEN

template <class ScalarImp>
class EigenRowMajorSparseMatrix
{
  static_assert(AlwaysFalse<ScalarImp>::value, "You are missing Eigen!");
};

template <class ScalarImp>
class EigenColMajorSparseMatrix
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
  using BaseType = LA::internal::MatrixAbstractionBase<LA::EigenRowMajorSparseMatrix<T>>;

  static constexpr Common::StorageLayout storage_layout = Common::StorageLayout::csr;

  template <size_t rows = BaseType::static_rows, size_t cols = BaseType::static_cols, class FieldType = T>
  using MatrixTypeTemplate = LA::EigenRowMajorSparseMatrix<FieldType>;
};

template <class T>
struct MatrixAbstraction<LA::EigenColMajorSparseMatrix<T>>
  : public LA::internal::MatrixAbstractionBase<LA::EigenColMajorSparseMatrix<T>>
{
  using BaseType = LA::internal::MatrixAbstractionBase<LA::EigenColMajorSparseMatrix<T>>;

  static constexpr Common::StorageLayout storage_layout = Common::StorageLayout::csc;

  template <size_t rows = BaseType::static_rows, size_t cols = BaseType::static_cols, class FieldType = T>
  using MatrixTypeTemplate = LA::EigenColMajorSparseMatrix<FieldType>;
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
