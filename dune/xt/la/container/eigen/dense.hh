// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// The copyright lies with the authors of this file (see below).
// License: Dual licensed as  BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Barbara Verfürth (2015)
//   Felix Schindler  (2014 - 2016)
//   Rene Milk        (2014 - 2015)
//   Tobias Leibner   (2014)

#ifndef DUNE_XT_LA_CONTAINER_EIGEN_DENSE_HH
#define DUNE_XT_LA_CONTAINER_EIGEN_DENSE_HH

#include <memory>
#include <type_traits>
#include <vector>
#include <initializer_list>
#include <complex>
#include <mutex>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/xt/common/disable_warnings.hh>
#if HAVE_EIGEN
#include <Eigen/Core>
#endif
#include <dune/xt/common/reenable_warnings.hh>

#include <dune/common/typetraits.hh>
#include <dune/common/densematrix.hh>
#include <dune/common/ftraits.hh>
#include <dune/common/unused.hh>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/crtp.hh>

#include "dune/xt/la/container/interfaces.hh"
#include "dune/xt/la/container/pattern.hh"

#include "base.hh"

namespace Dune {
namespace XT {
namespace LA {

// forwards
template <class ScalarImp>
class EigenDenseVector;

template <class T>
class EigenMappedDenseVector;

template <class ScalarImp>
class EigenDenseMatrix;

#if HAVE_EIGEN

namespace internal {

/**
 *  \brief Traits for EigenDenseVector.
 */
template <class ScalarImp = double>
class EigenDenseVectorTraits
{
public:
  typedef typename Dune::FieldTraits<ScalarImp>::field_type ScalarType;
  typedef typename Dune::FieldTraits<ScalarImp>::real_type RealType;
  typedef ScalarType DataType;
  typedef EigenDenseVector<ScalarType> derived_type;
  typedef typename ::Eigen::Matrix<ScalarType, ::Eigen::Dynamic, 1> BackendType;
  static const constexpr Backends backend_type = Backends::eigen_dense;
  static const constexpr Backends dense_matrix_type = Backends::eigen_dense;
  static const constexpr Backends sparse_matrix_type = Backends::eigen_sparse;
}; // class EigenDenseVectorTraits

/**
 *  \brief Traits for EigenMappedDenseVector.
 */
template <class ScalarImp = double>
class EigenMappedDenseVectorTraits
{
  typedef typename ::Eigen::Matrix<ScalarImp, ::Eigen::Dynamic, 1> PlainBackendType;

public:
  typedef typename Dune::FieldTraits<ScalarImp>::field_type ScalarType;
  typedef typename Dune::FieldTraits<ScalarImp>::real_type RealType;
  typedef EigenMappedDenseVector<ScalarType> derived_type;
  typedef Eigen::Map<PlainBackendType> BackendType;
  static const constexpr Backends backend_type = Backends::eigen_dense;
  static const constexpr Backends dense_matrix_type = Backends::eigen_dense;
  static const constexpr Backends sparse_matrix_type = Backends::eigen_sparse;
}; // class EigenMappedDenseVectorTraits

/**
 *  \brief Traits for EigenDenseMatrix.
 */
template <class ScalarImp = double>
class EigenDenseMatrixTraits
{
public:
  typedef typename Dune::FieldTraits<ScalarImp>::field_type ScalarType;
  typedef typename Dune::FieldTraits<ScalarImp>::real_type RealType;
  typedef ScalarType DataType;
  typedef EigenDenseMatrix<ScalarType> derived_type;
  typedef typename ::Eigen::Matrix<ScalarType, ::Eigen::Dynamic, ::Eigen::Dynamic> BackendType;
  static const constexpr Backends backend_type = Backends::eigen_dense;
  static const constexpr Backends vector_type = Backends::eigen_dense;
}; // class EigenDenseMatrixTraits

} // namespace internal

/**
 *  \brief A dense vector implementation of VectorInterface using the eigen backend.
 */
template <class ScalarImp = double>
class EigenDenseVector : public EigenBaseVector<internal::EigenDenseVectorTraits<ScalarImp>, ScalarImp>,
                         public ProvidesDataAccess<internal::EigenDenseVectorTraits<ScalarImp>>
{
  typedef EigenDenseVector<ScalarImp> ThisType;
  typedef VectorInterface<internal::EigenDenseVectorTraits<ScalarImp>, ScalarImp> VectorInterfaceType;
  typedef EigenBaseVector<internal::EigenDenseVectorTraits<ScalarImp>, ScalarImp> BaseType;

public:
  typedef internal::EigenDenseVectorTraits<ScalarImp> Traits;
  typedef typename Traits::ScalarType ScalarType;
  typedef typename Traits::RealType RealType;
  typedef typename Traits::DataType DataType;
  typedef typename Traits::BackendType BackendType;

private:
  typedef typename BackendType::Index EIGEN_size_t;

public:
  EigenDenseVector(const ThisType& other)
    : BaseType()
  {
    backend_ = other.backend_;
  }

  explicit EigenDenseVector(const size_t ss = 0, const ScalarType value = ScalarType(0))
  {
    backend_ = std::make_shared<BackendType>(ss);
    backend_->setOnes();
    backend_->operator*=(value);
  }

  explicit EigenDenseVector(const std::vector<ScalarType>& other)
  {
    backend_ = std::make_shared<BackendType>(internal::boost_numeric_cast<EIGEN_size_t>(other.size()));
    for (size_t ii = 0; ii < other.size(); ++ii)
      backend_->operator[](ii) = other[ii];
  }

  explicit EigenDenseVector(const std::initializer_list<ScalarType>& other)
  {
    backend_ = std::make_shared<BackendType>(internal::boost_numeric_cast<EIGEN_size_t>(other.size()));
    size_t ii = 0;
    for (auto element : other)
      backend_->operator[](ii++) = element;
  }

  explicit EigenDenseVector(const BackendType& other,
                            const bool /*prune*/ = false,
                            const ScalarType /*eps*/ = Common::FloatCmp::DefaultEpsilon<ScalarType>::value())
  {
    backend_ = std::make_shared<BackendType>(other);
  }

  /**
   *  \note Takes ownership of backend_ptr in the sense that you must not delete it afterwards!
   */
  explicit EigenDenseVector(BackendType* backend_ptr)
  {
    backend_ = std::shared_ptr<BackendType>(backend_ptr);
  }

  explicit EigenDenseVector(std::shared_ptr<BackendType> backend_ptr)
  {
    backend_ = backend_ptr;
  }

  using BaseType::operator=;

  /**
   *  \note Does a deep copy.
   */
  ThisType& operator=(const BackendType& other)
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(this->mutex_);
    backend_ = std::make_shared<BackendType>(other);
    return *this;
  } // ... operator=(...)

  using VectorInterfaceType::add;
  using VectorInterfaceType::sub;
  using BaseType::backend;

  /// \name Required by ProvidesDataAccess.
  /// \{

  /** \attention I have no idea if this makes sense for complex! **/
  DataType* data()
  {
    return backend().data();
  }

  size_t data_size() const
  {
    return this->size();
  }

  /// \}

private:
  using BaseType::backend_;

protected:
  inline void ensure_uniqueness()
  {
    if (!backend_.unique()) {
      std::lock_guard<std::mutex> DUNE_UNUSED(lock)(this->mutex_);
      backend_ = std::make_shared<BackendType>(*(backend_));
    }
  } // ... ensure_uniqueness(...)

private:
  friend class EigenBaseVector<internal::EigenDenseVectorTraits<ScalarType>, ScalarType>;
}; // class EigenDenseVector

/**
 *  \brief  A dense vector implementation of VectorInterface using the eigen backend which wrappes a raw array.
 */
template <class ScalarImp = double>
class EigenMappedDenseVector : public EigenBaseVector<internal::EigenMappedDenseVectorTraits<ScalarImp>, ScalarImp>
{
  typedef EigenMappedDenseVector<ScalarImp> ThisType;
  typedef VectorInterface<internal::EigenMappedDenseVectorTraits<ScalarImp>, ScalarImp> VectorInterfaceType;
  typedef EigenBaseVector<internal::EigenMappedDenseVectorTraits<ScalarImp>, ScalarImp> BaseType;
  static_assert(std::is_same<ScalarImp, double>::value, "Undefined behaviour for non-double data!");

public:
  typedef internal::EigenMappedDenseVectorTraits<ScalarImp> Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::ScalarType ScalarType;
  typedef typename Traits::RealType RealType;

private:
  typedef typename BackendType::Index EIGEN_size_t;

public:
  /**
   *  \brief  This is the constructor of interest which wrappes a raw array.
   */
  EigenMappedDenseVector(ScalarType* data, size_t data_size)
  {
    backend_ = std::make_shared<BackendType>(data, internal::boost_numeric_cast<EIGEN_size_t>(data_size));
  }

  /**
   *  \brief  This constructor allows to create an instance of this type just like any other vector.
   */
  explicit EigenMappedDenseVector(const size_t ss = 0, const ScalarType value = ScalarType(0))
  {
    backend_ = std::make_shared<BackendType>(new ScalarType[ss], internal::boost_numeric_cast<EIGEN_size_t>(ss));
    backend_->setOnes();
    backend_->operator*=(value);
  }

  explicit EigenMappedDenseVector(const std::vector<ScalarType>& other)
  {
    backend_ = std::make_shared<BackendType>(new ScalarType[other.size()],
                                             internal::boost_numeric_cast<EIGEN_size_t>(other.size()));
    for (size_t ii = 0; ii < other.size(); ++ii)
      backend_->operator[](ii) = other[ii];
  }

  explicit EigenMappedDenseVector(const std::initializer_list<ScalarType>& other)
  {
    backend_ = std::make_shared<BackendType>(new ScalarType[other.size()],
                                             internal::boost_numeric_cast<EIGEN_size_t>(other.size()));
    size_t ii = 0;
    for (auto element : other)
      backend_->operator[](ii++) = element;
  }

  /**
   *  \brief  This constructor does not do a deep copy.
   */
  EigenMappedDenseVector(const ThisType& other)
   : BaseType()
  {
    backend_ = other.backend_;
  }

  /**
   * \brief This constructor does a deep copy.
   */
  explicit EigenMappedDenseVector(const BackendType& other,
                                  const bool /*prune*/ = false,
                                  const ScalarType /*eps*/ = Common::FloatCmp::DefaultEpsilon<ScalarType>::value())
  {
    backend_ = std::make_shared<BackendType>(new ScalarType[other.size()],
                                             internal::boost_numeric_cast<EIGEN_size_t>(other.size()));
    backend_->operator=(other);
  }

  /**
   *  \note Takes ownership of backend_ptr in the sense that you must not delete it afterwards!
   */
  explicit EigenMappedDenseVector(BackendType* backend_ptr)
  {
    backend_ = std::shared_ptr<BackendType>(backend_ptr);
  }

  explicit EigenMappedDenseVector(std::shared_ptr<BackendType> backend_ptr)
  {
    backend_ = backend_ptr;
  }

  using BaseType::operator=;

  /**
   * \brief does a deep copy;
   */
  ThisType& operator=(const BackendType& other)
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(this->mutex_);
    backend_ = std::make_shared<BackendType>(new ScalarType[other.size()], other.size());
    backend_->operator=(other);
    return *this;
  }

  using VectorInterfaceType::add;
  using VectorInterfaceType::sub;
  using BaseType::backend;

private:
  using BaseType::backend_;

protected:
  inline void ensure_uniqueness()
  {
    if (!backend_.unique()) {
      std::lock_guard<std::mutex> DUNE_UNUSED(lock)(this->mutex_);
      auto new_backend = std::make_shared<BackendType>(new ScalarType[backend_->size()], backend_->size());
      new_backend->operator=(*(backend_));
      backend_ = new_backend;
    }
  } // ... ensure_uniqueness(...)

private:
  friend class EigenBaseVector<internal::EigenMappedDenseVectorTraits<ScalarType>, ScalarType>;
}; // class EigenMappedDenseVector

/**
 *  \brief  A dense matrix implementation of MatrixInterface using the eigen backend.
 */
template <class ScalarImp = double>
class EigenDenseMatrix : public MatrixInterface<internal::EigenDenseMatrixTraits<ScalarImp>, ScalarImp>,
                         public ProvidesBackend<internal::EigenDenseMatrixTraits<ScalarImp>>,
                         public ProvidesDataAccess<internal::EigenDenseMatrixTraits<ScalarImp>>
{
  typedef EigenDenseMatrix<ScalarImp> ThisType;
  typedef MatrixInterface<internal::EigenDenseMatrixTraits<ScalarImp>, ScalarImp> MatrixInterfaceType;

public:
  typedef internal::EigenDenseMatrixTraits<ScalarImp> Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::ScalarType ScalarType;
  typedef typename Traits::RealType RealType;
  typedef typename Traits::DataType DataType;

private:
  typedef typename BackendType::Index EIGEN_size_t;

public:
  explicit EigenDenseMatrix(const size_t rr = 0, const size_t cc = 0, const ScalarType value = ScalarType(0))
    : backend_(new BackendType(rr, cc))
  {
    this->backend_->setOnes();
    this->backend_->operator*=(value);
  }

  /// This constructors ignores the given pattern and initializes the matrix with 0.
  EigenDenseMatrix(const size_t rr, const size_t cc, const SparsityPatternDefault& /*pattern*/)
    : backend_(new BackendType(internal::boost_numeric_cast<EIGEN_size_t>(rr),
                               internal::boost_numeric_cast<EIGEN_size_t>(cc)))
  {
    backend_->setZero();
  }

  EigenDenseMatrix(const ThisType& other)
    : backend_(other.backend_)
  {
  }

  /**
   * \note If prune == true, this implementation is not optimal!
   */
  explicit EigenDenseMatrix(const BackendType& other,
                            const bool prune = false,
                            const typename Common::FloatCmp::DefaultEpsilon<ScalarType>::Type eps =
                                Common::FloatCmp::DefaultEpsilon<ScalarType>::value())
  {
    if (prune)
      backend_ = ThisType(other).pruned(eps).backend_;
    else
      backend_ = std::make_shared<BackendType>(other);
  }

  template <class M>
  EigenDenseMatrix(const MatrixInterface<M, ScalarType>& other)
    : backend_(new BackendType(other.rows(), other.cols()))
  {
    for (size_t ii = 0; ii < other.rows(); ++ii)
      for (size_t jj = 0; jj < other.cols(); ++jj)
        set_entry(ii, jj, other.get_entry(ii, jj));
  }

  template <class T>
  EigenDenseMatrix(const DenseMatrix<T>& other)
    : backend_(new BackendType(other.rows(), other.cols()))
  {
    for (size_t ii = 0; ii < other.rows(); ++ii)
      for (size_t jj = 0; jj < other.cols(); ++jj)
        set_entry(ii, jj, other[ii][jj]);
  }

  /**
   *  \note Takes ownership of backend_ptr in the sense that you must not delete it afterwards!
   */
  explicit EigenDenseMatrix(BackendType* backend_ptr)
    : backend_(backend_ptr)
  {
  }

  explicit EigenDenseMatrix(std::shared_ptr<BackendType> backend_ptr)
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
  /// \name Required by the ProvidesDataAccess interface.
  /// \{

  DataType* data()
  {
    return backend_->data();
  }

  size_t data_size() const
  {
    return rows() * cols();
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
    auto& backend_ref = backend();
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    backend_ref *= alpha;
  }

  void axpy(const ScalarType& alpha, const ThisType& xx)
  {
    auto& backend_ref = backend();
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
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
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    yy.backend().transpose() = backend() * xx.backend();
  }

  void add_to_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    auto& backend_ref = backend();
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    assert(ii < rows());
    assert(jj < cols());
    backend_ref(ii, jj) += value;
  } // ... add_to_entry(...)

  void set_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    auto& backend_ref = backend();
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    assert(ii < rows());
    assert(jj < cols());
    backend_ref(ii, jj) = value;
  } // ... set_entry(...)

  ScalarType get_entry(const size_t ii, const size_t jj) const
  {
    assert(ii < rows());
    assert(jj < cols());
    return backend()(ii, jj);
  } // ... get_entry(...)

  void clear_row(const size_t ii)
  {
    auto& backend_ref = backend();
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    if (ii >= rows())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given ii (" << ii << ") is larger than the rows of this (" << rows() << ")!");
    for (size_t jj = 0; jj < cols(); ++jj)
      backend_ref(ii, jj) = ScalarType(0);
  } // ... clear_row(...)

  void clear_col(const size_t jj)
  {
    auto& backend_ref = backend();
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    if (jj >= cols())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given jj (" << jj << ") is larger than the cols of this (" << cols() << ")!");
    for (size_t ii = 0; ii < rows(); ++ii)
      backend_ref(ii, jj) = ScalarType(0);
  } // ... clear_col(...)

  void unit_row(const size_t ii)
  {
    auto& backend_ref = backend();
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    if (ii >= cols())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given ii (" << ii << ") is larger than the cols of this (" << cols() << ")!");
    if (ii >= rows())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given ii (" << ii << ") is larger than the rows of this (" << rows() << ")!");
    for (size_t jj = 0; jj < cols(); ++jj)
      backend_ref(ii, jj) = ScalarType(0);
    backend_ref(ii, ii) = ScalarType(1);
  } // ... unit_row(...)

  void unit_col(const size_t jj)
  {
    auto& backend_ref = backend();
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    if (jj >= cols())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given jj (" << jj << ") is larger than the cols of this (" << cols() << ")!");
    if (jj >= rows())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given jj (" << jj << ") is larger than the rows of this (" << rows() << ")!");
    for (size_t ii = 0; ii < rows(); ++ii)
      backend_ref(ii, jj) = ScalarType(0);
    backend_ref(jj, jj) = ScalarType(1);
  } // ... unit_col(...)

  bool valid() const
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
    for (size_t ii = 0; ii < rows(); ++ii) {
      for (size_t jj = 0; jj < cols(); ++jj) {
        const auto& entry = backend()(ii, jj);
        if (Common::isnan(entry) || Common::isinf(entry))
          return false;
      }
    }
    return true;
  } // ... valid(...)

  /**
   * \}
   */

protected:
  inline void ensure_uniqueness()
  {
    if (!backend_.unique()) {
      std::lock_guard<std::mutex> DUNE_UNUSED(lock)(mutex_);
      backend_ = std::make_shared<BackendType>(*backend_);
    }
  } // ... ensure_uniqueness(...)

private:
  std::shared_ptr<BackendType> backend_;
  mutable std::mutex mutex_;
}; // class EigenDenseMatrix

#else // HAVE_EIGEN

template <class ScalarImp>
class EigenDenseVector
{
  static_assert(Dune::AlwaysFalse<ScalarImp>::value, "You are missing Eigen!");
};

template <class ScalarImp>
class EigenMappedDenseVector
{
  static_assert(Dune::AlwaysFalse<ScalarImp>::value, "You are missing Eigen!");
};

template <class ScalarImp>
class EigenDenseMatrix
{
  static_assert(Dune::AlwaysFalse<ScalarImp>::value, "You are missing Eigen!");
};

#endif // HAVE_EIGEN

} // namespace LA
namespace Common {

#if HAVE_EIGEN

template <class T>
struct VectorAbstraction<LA::EigenDenseVector<T>> : public LA::internal::VectorAbstractionBase<LA::EigenDenseVector<T>>
{
};

template <class T>
struct VectorAbstraction<LA::EigenMappedDenseVector<T>>
    : public LA::internal::VectorAbstractionBase<LA::EigenMappedDenseVector<T>>
{
};

template <class T>
struct MatrixAbstraction<LA::EigenDenseMatrix<T>> : public LA::internal::MatrixAbstractionBase<LA::EigenDenseMatrix<T>>
{
};

#endif // HAVE_EIGEN

} // namespace Common
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_CONTAINER_EIGEN_DENSE_HH
