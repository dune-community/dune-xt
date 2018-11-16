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
//   Tobias Leibner   (2014, 2016 - 2018)

#ifndef DUNE_XT_LA_CONTAINER_EIGEN_DENSE_HH
#define DUNE_XT_LA_CONTAINER_EIGEN_DENSE_HH

#include <memory>
#include <type_traits>
#include <vector>
#include <initializer_list>
#include <complex>
#include <mutex>

#include <boost/numeric/conversion/cast.hpp>

#if HAVE_EIGEN
#include <dune/xt/common/disable_warnings.hh>
#include <Eigen/Core>
#include <dune/xt/common/reenable_warnings.hh>
#endif

#include <dune/common/typetraits.hh>
#include <dune/common/densematrix.hh>
#include <dune/common/ftraits.hh>
#include <dune/common/unused.hh>

#include <dune/xt/common/crtp.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/numeric_cast.hh>

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
class EigenDenseVectorTraits : public VectorTraitsBase<ScalarImp,
                                                       EigenDenseVector<ScalarImp>,
                                                       ::Eigen::Matrix<ScalarImp, ::Eigen::Dynamic, 1>,
                                                       Backends::eigen_dense,
                                                       Backends::eigen_dense,
                                                       Backends::eigen_sparse>
{
};


/**
 *  \brief Traits for EigenMappedDenseVector.
 */
template <class ScalarImp = double>
class EigenMappedDenseVectorTraits
    : public VectorTraitsBase<ScalarImp,
                              EigenMappedDenseVector<ScalarImp>,
                              Eigen::Map<::Eigen::Matrix<ScalarImp, ::Eigen::Dynamic, 1>>,
                              Backends::eigen_dense,
                              Backends::eigen_dense,
                              Backends::eigen_sparse>
{
};


/**
 *  \brief Traits for EigenDenseMatrix.
 */
template <class ScalarImp = double>
class EigenDenseMatrixTraits : public MatrixTraitsBase<ScalarImp,
                                                       EigenDenseMatrix<ScalarImp>,
                                                       ::Eigen::Matrix<ScalarImp, ::Eigen::Dynamic, ::Eigen::Dynamic>,
                                                       Backends::eigen_dense,
                                                       Backends::eigen_dense>
{
};


} // namespace internal


/**
 *  \brief A dense vector implementation of VectorInterface using the eigen backend.
 */
template <class ScalarImp = double>
class EigenDenseVector : public EigenBaseVector<internal::EigenDenseVectorTraits<ScalarImp>, ScalarImp>,
                         public ProvidesDataAccess<internal::EigenDenseVectorTraits<ScalarImp>>
{
  using ThisType = EigenDenseVector;
  using InterfaceType = VectorInterface<internal::EigenDenseVectorTraits<ScalarImp>, ScalarImp>;
  using BaseType = EigenBaseVector<internal::EigenDenseVectorTraits<ScalarImp>, ScalarImp>;

public:
  using typename InterfaceType::RealType;
  using typename InterfaceType::ScalarType;
  using Traits = typename InterfaceType::Traits;
  using typename ProvidesBackend<Traits>::BackendType;
  using typename ProvidesDataAccess<Traits>::DataType;
  //! derived_type is ambiguous due to multiple parent classes otherwise
  using derived_type = typename InterfaceType::derived_type;

private:
  using EIGEN_size_t = typename BackendType::Index;

public:
  EigenDenseVector(const ThisType& other) = default;

  explicit EigenDenseVector(const size_t ss = 0, const ScalarType value = ScalarType(0), const size_t num_mutexes = 1)
    : BaseType(num_mutexes)
  {
    backend_ = std::make_shared<BackendType>(ss);
    backend_->setOnes();
    backend_->operator*=(value);
  }

  explicit EigenDenseVector(const std::vector<ScalarType>& other, const size_t num_mutexes = 1)
    : BaseType(num_mutexes)
  {
    backend_ = std::make_shared<BackendType>(Common::numeric_cast<EIGEN_size_t>(other.size()));
    for (size_t ii = 0; ii < other.size(); ++ii)
      backend_->operator[](ii) = other[ii];
  }

  explicit EigenDenseVector(const std::initializer_list<ScalarType>& other, const size_t num_mutexes = 1)
    : BaseType(num_mutexes)
  {
    backend_ = std::make_shared<BackendType>(Common::numeric_cast<EIGEN_size_t>(other.size()));
    size_t ii = 0;
    for (auto element : other)
      backend_->operator[](ii++) = element;
  }

  explicit EigenDenseVector(const BackendType& other,
                            const bool /*prune*/ = false,
                            const ScalarType /*eps*/ = Common::FloatCmp::DefaultEpsilon<ScalarType>::value(),
                            const size_t num_mutexes = 1)
    : BaseType(num_mutexes)
  {
    backend_ = std::make_shared<BackendType>(other);
  }

  /**
   *  \note Takes ownership of backend_ptr in the sense that you must not delete it afterwards!
   */
  explicit EigenDenseVector(BackendType* backend_ptr, const size_t num_mutexes = 1)
    : BaseType(num_mutexes)
  {
    backend_ = std::shared_ptr<BackendType>(backend_ptr);
  }

  explicit EigenDenseVector(std::shared_ptr<BackendType> backend_ptr, const size_t num_mutexes = 1)
    : BaseType(num_mutexes)
  {
    backend_ = backend_ptr;
  }

  using BaseType::operator=;

  ThisType& operator=(const BackendType& other)
  {
    backend_ = std::make_shared<BackendType>(other);
    return *this;
  } // ... operator=(...)

  using InterfaceType::add;
  using InterfaceType::sub;
  using InterfaceType::operator+;
  using InterfaceType::operator-;
  using InterfaceType::operator*;
  using BaseType::backend;

  /// \name Required by VectorInterface
  /// \{

  inline void resize(const size_t new_size)
  {
    // Eigens resize loses all values, so we have to store them before...
    auto stored_size = new_size > this->size() ? this->size() : new_size;
    std::vector<ScalarType> old_values(stored_size);
    for (size_t ii = 0; ii < stored_size; ++ii)
      old_values[ii] = this->get_entry(ii);
    backend_->resize(new_size, 1);
    for (size_t ii = 0; ii < stored_size; ++ii)
      this->set_entry(ii, old_values[ii]);
  }

  /// \}
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

  friend class VectorInterface<internal::EigenDenseVectorTraits<ScalarType>, ScalarType>;
  friend class EigenBaseVector<internal::EigenDenseVectorTraits<ScalarType>, ScalarType>;
}; // class EigenDenseVector

/**
 *  \brief  A dense vector implementation of VectorInterface using the eigen backend which wrappes a raw array.
 */
template <class ScalarImp = double>
class EigenMappedDenseVector : public EigenBaseVector<internal::EigenMappedDenseVectorTraits<ScalarImp>, ScalarImp>
{
  using ThisType = EigenMappedDenseVector;
  using InterfaceType = VectorInterface<internal::EigenMappedDenseVectorTraits<ScalarImp>, ScalarImp>;
  using BaseType = EigenBaseVector<internal::EigenMappedDenseVectorTraits<ScalarImp>, ScalarImp>;
  static_assert(std::is_same<ScalarImp, double>::value, "Undefined behaviour for non-double data!");

public:
  using typename InterfaceType::RealType;
  using typename InterfaceType::ScalarType;
  using Traits = typename InterfaceType::Traits;
  using BackendType = typename Traits::BackendType;
  // needed to fix gcc compilation error due to ambiguous lookup of derived type
  using derived_type = typename Traits::derived_type;

private:
  using EIGEN_size_t = typename BackendType::Index;

public:
  /**
   *  \brief  This is the constructor of interest which wrappes a raw array.
   */
  EigenMappedDenseVector(ScalarType* data, size_t data_size, const size_t num_mutexes = 1)
    : BaseType(num_mutexes)
  {
    backend_ = std::make_shared<BackendType>(data, Common::numeric_cast<EIGEN_size_t>(data_size));
  }

  /**
   *  \brief  This constructor allows to create an instance of this type just like any other vector.
   */
  explicit EigenMappedDenseVector(const size_t ss = 0,
                                  const ScalarType value = ScalarType(0),
                                  const size_t num_mutexes = 1)
    : BaseType(num_mutexes)
  {
    backend_ = std::make_shared<BackendType>(new ScalarType[ss], Common::numeric_cast<EIGEN_size_t>(ss));
    backend_->setOnes();
    backend_->operator*=(value);
  }

  explicit EigenMappedDenseVector(const std::vector<ScalarType>& other, const size_t num_mutexes = 1)
    : BaseType(num_mutexes)
  {
    backend_ =
        std::make_shared<BackendType>(new ScalarType[other.size()], Common::numeric_cast<EIGEN_size_t>(other.size()));
    for (size_t ii = 0; ii < other.size(); ++ii)
      (*backend_)[ii] = other[ii];
  }

  explicit EigenMappedDenseVector(const std::initializer_list<ScalarType>& other, const size_t num_mutexes = 1)
    : BaseType(num_mutexes)
  {
    backend_ =
        std::make_shared<BackendType>(new ScalarType[other.size()], Common::numeric_cast<EIGEN_size_t>(other.size()));
    size_t ii = 0;
    for (auto element : other)
      (*backend_)[ii++] = element;
  }

  EigenMappedDenseVector(const ThisType& other)
    : BaseType(other)
  {
    backend_ =
        std::make_shared<BackendType>(new ScalarType[other.size()], Common::numeric_cast<EIGEN_size_t>(other.size()));
    backend_->operator=(other.backend());
  }

  explicit EigenMappedDenseVector(const BackendType& other,
                                  const bool /*prune*/ = false,
                                  const ScalarType /*eps*/ = Common::FloatCmp::DefaultEpsilon<ScalarType>::value(),
                                  const size_t num_mutexes = 1)
    : BaseType(num_mutexes)
  {
    backend_ =
        std::make_shared<BackendType>(new ScalarType[other.size()], Common::numeric_cast<EIGEN_size_t>(other.size()));
    backend_->operator=(other);
  }

  /**
   *  \note Takes ownership of backend_ptr in the sense that you must not delete it afterwards!
   */
  explicit EigenMappedDenseVector(BackendType* backend_ptr, const size_t num_mutexes = 1)
    : BaseType(num_mutexes)
  {
    backend_ = std::shared_ptr<BackendType>(backend_ptr);
  }

  explicit EigenMappedDenseVector(std::shared_ptr<BackendType> backend_ptr, const size_t num_mutexes = 1)
    : BaseType(num_mutexes)
  {
    backend_ = backend_ptr;
  }

  using BaseType::operator=;

  ThisType& operator=(const ThisType& other)
  {
    BaseType::operator=(other);
    return *this;
  }

  ThisType& operator=(const BackendType& other)
  {
    backend_ = std::make_shared<BackendType>(new ScalarType[other.size()], other.size());
    backend_->operator=(other);
    return *this;
  }

  inline void resize(const size_t new_size)
  {
    if (new_size != this->size())
      DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong, "This does not make sense for mapped memory!");
  }

  using InterfaceType::add;
  using InterfaceType::sub;
  using InterfaceType::operator+;
  using InterfaceType::operator-;
  using InterfaceType::operator*;
  using BaseType::backend;

private:
  using BaseType::backend_;

  friend class VectorInterface<internal::EigenMappedDenseVectorTraits<ScalarType>, ScalarType>;
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

  using ThisType = EigenDenseMatrix;
  using InterfaceType = MatrixInterface<internal::EigenDenseMatrixTraits<ScalarImp>, ScalarImp>;

public:
  using typename InterfaceType::RealType;
  using typename InterfaceType::ScalarType;
  using Traits = typename InterfaceType::Traits;
  using typename ProvidesBackend<Traits>::BackendType;
  using typename ProvidesDataAccess<Traits>::DataType;
  //! derived_type is ambiguous due to multiple parent classes otherwise
  using derived_type = typename InterfaceType::derived_type;

private:
  using MutexesType = typename Traits::MutexesType;
  using EIGEN_size_t = typename BackendType::Index;

public:
  explicit EigenDenseMatrix(const size_t rr = 0,
                            const size_t cc = 0,
                            const ScalarType value = ScalarType(0),
                            const size_t num_mutexes = 1)
    : backend_(new BackendType(rr, cc))
    , mutexes_(std::make_unique<MutexesType>(num_mutexes))
  {
    this->backend_->setOnes();
    this->backend_->operator*=(value);
  }

  /// This constructors ignores the given pattern and initializes the matrix with 0.
  EigenDenseMatrix(const size_t rr,
                   const size_t cc,
                   const SparsityPatternDefault& /*pattern*/,
                   const size_t num_mutexes = 1)
    : backend_(new BackendType(Common::numeric_cast<EIGEN_size_t>(rr), Common::numeric_cast<EIGEN_size_t>(cc)))
    , mutexes_(std::make_unique<MutexesType>(num_mutexes))
  {
    backend_->setZero();
  }

  EigenDenseMatrix(const ThisType& other)
    : backend_(std::make_shared<BackendType>(*other.backend_))
    , mutexes_(std::make_unique<MutexesType>(other.mutexes_->size()))
  {
  }

  /**
   * \note If prune == true, this implementation is not optimal!
   */
  explicit EigenDenseMatrix(const BackendType& other,
                            const bool prune = false,
                            const typename Common::FloatCmp::DefaultEpsilon<ScalarType>::Type eps =
                                Common::FloatCmp::DefaultEpsilon<ScalarType>::value(),
                            const size_t num_mutexes = 1)
    : mutexes_(std::make_unique<MutexesType>(num_mutexes))
  {
    if (prune)
      backend_ = ThisType(other).pruned(eps).backend_;
    else
      backend_ = std::make_shared<BackendType>(other);
  }

  template <class M>
  EigenDenseMatrix(const MatrixInterface<M, ScalarType>& other, const size_t num_mutexes = 1)
    : backend_(new BackendType(other.rows(), other.cols()))
    , mutexes_(std::make_unique<MutexesType>(num_mutexes))
  {
    for (size_t ii = 0; ii < other.rows(); ++ii)
      for (size_t jj = 0; jj < other.cols(); ++jj)
        set_entry(ii, jj, other.get_entry(ii, jj));
  }

  template <class T>
  EigenDenseMatrix(const DenseMatrix<T>& other, const size_t num_mutexes = 1)
    : backend_(new BackendType(other.rows(), other.cols()))
    , mutexes_(std::make_unique<MutexesType>(num_mutexes))
  {
    for (size_t ii = 0; ii < other.rows(); ++ii)
      for (size_t jj = 0; jj < other.cols(); ++jj)
        set_entry(ii, jj, other[ii][jj]);
  }

  /**
   *  \note Takes ownership of backend_ptr in the sense that you must not delete it afterwards!
   */
  explicit EigenDenseMatrix(BackendType* backend_ptr, const size_t num_mutexes = 1)
    : backend_(backend_ptr)
    , mutexes_(std::make_unique<MutexesType>(num_mutexes))
  {
  }

  explicit EigenDenseMatrix(std::shared_ptr<BackendType> backend_ptr, const size_t num_mutexes = 1)
    : backend_(backend_ptr)
    , mutexes_(std::make_unique<MutexesType>(num_mutexes))
  {
  }

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      *backend_ = *other.backend_;
      mutexes_ = std::make_unique<MutexesType>(other.mutexes_->size());
    }
    return *this;
  }

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
    return ThisType(*backend_);
  }

  void scal(const ScalarType& alpha)
  {
    const internal::VectorLockGuard DUNE_UNUSED(guard)(*mutexes_);
    backend() *= alpha;
  }

  void axpy(const ScalarType& alpha, const ThisType& xx)
  {
    const internal::VectorLockGuard DUNE_UNUSED(guard)(*mutexes_);
    if (!has_equal_shape(xx))
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The shape of xx (" << xx.rows() << "x" << xx.cols() << ") does not match the shape of this ("
                                     << rows()
                                     << "x"
                                     << cols()
                                     << ")!");
    backend() += alpha * xx.backend();
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
    internal::LockGuard DUNE_UNUSED(lock)(*mutexes_, jj, cols());
    assert(ii < rows());
    assert(jj < cols());
    backend()(ii, jj) += value;
  } // ... add_to_entry(...)

  void set_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    assert(ii < rows());
    assert(jj < cols());
    backend()(ii, jj) = value;
  } // ... set_entry(...)

  ScalarType get_entry(const size_t ii, const size_t jj) const
  {
    assert(ii < rows());
    assert(jj < cols());
    return backend()(ii, jj);
  } // ... get_entry(...)

  void clear_row(const size_t ii)
  {
    if (ii >= rows())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given ii (" << ii << ") is larger than the rows of this (" << rows() << ")!");
    for (size_t jj = 0; jj < cols(); ++jj)
      backend()(ii, jj) = ScalarType(0);
  } // ... clear_row(...)

  void clear_col(const size_t jj)
  {
    if (jj >= cols())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given jj (" << jj << ") is larger than the cols of this (" << cols() << ")!");
    for (size_t ii = 0; ii < rows(); ++ii)
      backend()(ii, jj) = ScalarType(0);
  } // ... clear_col(...)

  void unit_row(const size_t ii)
  {
    if (ii >= cols())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given ii (" << ii << ") is larger than the cols of this (" << cols() << ")!");
    if (ii >= rows())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given ii (" << ii << ") is larger than the rows of this (" << rows() << ")!");
    for (size_t jj = 0; jj < cols(); ++jj)
      backend()(ii, jj) = ScalarType(0);
    backend()(ii, ii) = ScalarType(1);
  } // ... unit_row(...)

  void unit_col(const size_t jj)
  {
    if (jj >= cols())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given jj (" << jj << ") is larger than the cols of this (" << cols() << ")!");
    if (jj >= rows())
      DUNE_THROW(Common::Exceptions::index_out_of_range,
                 "Given jj (" << jj << ") is larger than the rows of this (" << rows() << ")!");
    for (size_t ii = 0; ii < rows(); ++ii)
      backend()(ii, jj) = ScalarType(0);
    backend()(jj, jj) = ScalarType(1);
  } // ... unit_col(...)

  bool valid() const
  {
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
   * \{
   * \name Provided for improved performance.
   */

  using InterfaceType::operator-;

  ThisType operator-(const ThisType& other) const override final
  {
    return ThisType(this->backend() - other.backend());
  }

  using InterfaceType::operator+;

  ThisType operator+(const ThisType& other) const override final
  {
    return ThisType(this->backend() + other.backend());
  }

  using InterfaceType::operator*;

  ThisType operator*(const ThisType& other) const override final
  {
    return ThisType(this->backend() * other.backend());
  }

  /**
   * \}
   */
private:
  std::shared_ptr<BackendType> backend_;
  std::unique_ptr<std::vector<std::mutex>> mutexes_;
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
  using BaseType = LA::internal::MatrixAbstractionBase<LA::EigenDenseMatrix<T>>;

  static constexpr XT::Common::StorageLayout storage_layout = XT::Common::StorageLayout::dense_column_major;

  template <size_t rows = BaseType::static_rows, size_t cols = BaseType::static_cols, class FieldType = T>
  using MatrixTypeTemplate = LA::EigenDenseMatrix<FieldType>;

  static inline T* data(LA::EigenDenseMatrix<T>& mat)
  {
    return &(mat.backend()(0, 0));
  }

  static inline const T* data(const LA::EigenDenseMatrix<T>& mat)
  {
    return &(mat.backend()(0, 0));
  }
};

#endif // HAVE_EIGEN

} // namespace Common
} // namespace XT
} // namespace Dune


// begin: this is what we need for the lib
#if DUNE_XT_WITH_PYTHON_BINDINGS && HAVE_EIGEN


extern template class Dune::XT::LA::EigenDenseVector<double>;
extern template class Dune::XT::LA::EigenMappedDenseVector<double>;
extern template class Dune::XT::LA::EigenDenseMatrix<double>;
// extern template void Dune::XT::LA::EigenDenseMatrix<double>::mv(const Dune::XT::LA::EigenDenseVector<double>&,
//                                                                Dune::XT::LA::EigenDenseVector<double>&) const;
// extern template void Dune::XT::LA::EigenDenseMatrix<double>::mv(const Dune::XT::LA::EigenMappedDenseVector<double>&,
//                                                                Dune::XT::LA::EigenMappedDenseVector<double>&) const;
// extern template void Dune::XT::LA::EigenDenseMatrix<double>::mv(const Dune::XT::LA::EigenDenseVector<double>&,
//                                                                Dune::XT::LA::EigenMappedDenseVector<double>&) const;
// extern template void Dune::XT::LA::EigenDenseMatrix<double>::mv(const Dune::XT::LA::EigenMappedDenseVector<double>&,
//                                                                Dune::XT::LA::EigenDenseVector<double>&) const;


#endif // DUNE_XT_WITH_PYTHON_BINDINGS && HAVE_EIGEN
// end: this is what we need for the lib


#endif // DUNE_XT_LA_CONTAINER_EIGEN_DENSE_HH
