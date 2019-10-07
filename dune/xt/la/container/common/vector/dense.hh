// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)
//   Tobias Leibner  (2017 - 2019)

#ifndef DUNE_XT_LA_CONTAINER_COMMON_VECTOR_DENSE_HH
#define DUNE_XT_LA_CONTAINER_COMMON_VECTOR_DENSE_HH

#include <algorithm>
#include <cmath>
#include <complex>
#include <initializer_list>
#include <memory>
#include <mutex>
#include <numeric>
#include <type_traits>
#include <vector>

#include <dune/common/dynvector.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/ftraits.hh>
#include <dune/common/unused.hh>

#include <dune/xt/common/exceptions.hh>

#include <dune/xt/la/container/vector-interface.hh>

namespace Dune {
namespace XT {
namespace LA {


// forwards
template <class ScalarImp>
class CommonDenseVector;


namespace internal {


template <class ScalarImp>
struct CommonDenseVectorBackend
{
  using ScalarType = ScalarImp;
  using RealType = typename Dune::FieldTraits<ScalarImp>::real_type;
  using ThisType = CommonDenseVectorBackend;

  CommonDenseVectorBackend(const size_t ss = 0, const ScalarType& value = ScalarType())
    : size_(ss)
    , values_vector_(size_, value)
    , values_ptr_(values_vector_.data())
    , created_from_ptr_(false)
  {}

  CommonDenseVectorBackend(const size_t ss, ScalarType* values_ptr)
    : size_(ss)
    , values_ptr_(values_ptr)
    , created_from_ptr_(true)
  {}

  CommonDenseVectorBackend(const ThisType& other)
    : size_(other.size_)
    , values_vector_(other.size_)
    , values_ptr_(values_vector_.data())
    , created_from_ptr_(false)
  {
    std::copy_n(
        other.created_from_ptr_ ? other.values_ptr_ : other.values_vector_.data(), size_, values_vector_.begin());
  }

  CommonDenseVectorBackend(ThisType&& other)
    : size_(std::move(other.size_))
    , values_vector_(std::move(other.values_vector_))
    , values_ptr_(other.created_from_ptr_ ? other.values_ptr_ : values_vector_.data())
    , created_from_ptr_(other.created_from_ptr_)
  {}

  ThisType copy()
  {
    return ThisType(*this);
  }

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      size_ = other.size_;
      if (other.created_from_ptr_) {
        values_vector_ = std::vector<ScalarType>(size_);
        std::copy_n(other.values_ptr_, size_, values_vector_.begin());
      } else {
        values_vector_ = other.values_vector_;
      }
      values_ptr_ = values_vector_.data();
    }
    return *this;
  }

  ScalarType& operator[](const size_t ii)
  {
    return values_ptr_[ii];
  }

  const ScalarType& operator[](const size_t ii) const
  {
    return values_ptr_[ii];
  }

  void resize(const size_t ss)
  {
    size_ = ss;
    values_vector_.resize(size_);
    values_ptr_ = values_vector_.data();
  }

  size_t size() const
  {
    return size_;
  }

  ThisType& operator+=(const ThisType& other)
  {
    assert(other.size_ == size_);
    for (size_t ii = 0; ii < size_; ++ii)
      values_ptr_[ii] += other[ii];
    return *this;
  }

  ThisType& operator-=(const ThisType& other)
  {
    assert(other.size_ == size_);
    for (size_t ii = 0; ii < size_; ++ii)
      values_ptr_[ii] -= other[ii];
    return *this;
  }

  ThisType& operator*=(const ThisType& other)
  {
    assert(other.size_ == size_);
    for (size_t ii = 0; ii < size_; ++ii)
      values_ptr_[ii] *= other[ii];
    return *this;
  }

  ThisType& operator*=(const ScalarType& val)
  {
    for (size_t ii = 0; ii < size_; ++ii)
      values_ptr_[ii] *= val;
    return *this;
  }

  ThisType& operator/=(const ThisType& other)
  {
    assert(other.size_ == size_);
    for (size_t ii = 0; ii < size_; ++ii)
      values_ptr_[ii] /= other[ii];
    return *this;
  }

  ScalarType operator*(const ThisType& other) const
  {
    return std::inner_product(values_ptr_, values_ptr_ + size_, other.values_ptr_, ScalarType(0.));
  }

  RealType l1_norm() const
  {
    return std::accumulate(values_ptr_, values_ptr_ + size_, RealType(0.), [](const RealType& a, const ScalarType& b) {
      return a + std::abs(b);
    });
  }

  RealType l2_norm() const
  {
    return std::sqrt(
        std::accumulate(values_ptr_, values_ptr_ + size_, RealType(0.), [](const RealType& a, const ScalarType& b) {
          return a + std::pow(std::abs(b), 2);
        }));
  }

  RealType sup_norm() const
  {
    return std::abs(*std::max_element(values_ptr_, values_ptr_ + size_, [](const ScalarType& a, const ScalarType& b) {
      return std::abs(a) < std::abs(b);
    }));
  }

  size_t size_;
  std::vector<ScalarType> values_vector_;
  ScalarType* values_ptr_;
  bool created_from_ptr_;
};

/// Traits for CommonDenseVector
template <class ScalarImp>
class CommonDenseVectorTraits
  : public VectorTraitsBase<ScalarImp,
                            CommonDenseVector<ScalarImp>,
                            CommonDenseVectorBackend<ScalarImp>,
                            Backends::common_dense,
                            Backends::common_dense,
                            Backends::common_sparse>
{};


} // namespace internal


/**
 *  \brief A dense vector implementation of VectorInterface using the Dune::DynamicVector.
 */
template <class ScalarImp = double>
class CommonDenseVector
  : public VectorInterface<internal::CommonDenseVectorTraits<ScalarImp>, ScalarImp>
  , public ProvidesBackend<internal::CommonDenseVectorTraits<ScalarImp>>
  , public ProvidesDataAccess<internal::CommonDenseVectorTraits<ScalarImp>>
{
  using ThisType = CommonDenseVector;
  using InterfaceType = VectorInterface<internal::CommonDenseVectorTraits<ScalarImp>, ScalarImp>;

public:
  using typename InterfaceType::RealType;
  using typename InterfaceType::ScalarType;
  using Traits = typename InterfaceType::Traits;
  using typename ProvidesBackend<Traits>::BackendType;
  using typename ProvidesDataAccess<Traits>::DataType;
  // needed to fix gcc compilation error due to ambiguous lookup of derived type
  using derived_type = typename Traits::derived_type;

private:
  using MutexesType = typename Traits::MutexesType;

public:
  explicit CommonDenseVector(const size_t ss = 0, const ScalarType& value = ScalarType(), const size_t num_mutexes = 1)
    : backend_(new BackendType(ss, value))
    , mutexes_(std::make_unique<MutexesType>(num_mutexes))
  {}

  explicit CommonDenseVector(const std::vector<ScalarType>& other, const size_t num_mutexes = 1)
    : backend_(new BackendType(other.size()))
    , mutexes_(std::make_unique<MutexesType>(num_mutexes))
  {
    for (size_t ii = 0; ii < other.size(); ++ii)
      backend_->operator[](ii) = other[ii];
  }

  explicit CommonDenseVector(const std::initializer_list<ScalarType>& other, const size_t num_mutexes = 1)
    : backend_(new BackendType(other.size()))
    , mutexes_(std::make_unique<MutexesType>(num_mutexes))
  {
    size_t ii = 0;
    for (auto element : other) {
      backend_->operator[](ii) = element;
      ++ii;
    }
  } // CommonDenseVector(...)

  explicit CommonDenseVector(const size_t ss, ScalarType* values_ptr, const size_t num_mutexes = 1)
    : backend_(new BackendType(ss, values_ptr))
    , mutexes_(std::make_unique<MutexesType>(num_mutexes))
  {} // CommonDenseVector(...)

  explicit CommonDenseVector(const Dune::DynamicVector<ScalarType>& dynvec, const size_t num_mutexes = 1)
    : backend_(new BackendType(dynvec.size()))
    , mutexes_(std::make_unique<MutexesType>(num_mutexes))
  {
    for (size_t ii = 0; ii < size(); ++ii)
      (*backend_)[ii] = dynvec[ii];
  } // CommonDenseVector(...)

  CommonDenseVector(const ThisType& other)
    : backend_(std::make_shared<BackendType>(*other.backend_))
    , mutexes_(std::make_unique<MutexesType>(other.mutexes_->size()))
  {}

  explicit CommonDenseVector(const BackendType& other,
                             const bool /*prune*/ = false,
                             const ScalarType /*eps*/ = Common::FloatCmp::DefaultEpsilon<ScalarType>::value(),
                             const size_t num_mutexes = 1)
    : backend_(new BackendType(other))
    , mutexes_(std::make_unique<MutexesType>(num_mutexes))
  {}

  /**
   *  \note Takes ownership of backend_ptr in the sense that you must not delete it afterwards!
   */
  explicit CommonDenseVector(BackendType* backend_ptr, const size_t num_mutexes = 1)
    : backend_(backend_ptr)
    , mutexes_(std::make_unique<MutexesType>(num_mutexes))
  {}

  explicit CommonDenseVector(std::shared_ptr<BackendType> backend_ptr, const size_t num_mutexes = 1)
    : backend_(backend_ptr)
    , mutexes_(std::make_unique<MutexesType>(num_mutexes))
  {}

  using InterfaceType::operator=;

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      *backend_ = *other.backend_;
      mutexes_ = std::make_unique<MutexesType>(other.mutexes_->size());
    }
    return *this;
  }

  ThisType& operator=(const ScalarType& value)
  {
    for (auto& element : *this)
      element = value;
    return *this;
  } // ... operator=(...)

  /**
   *  \note Does a deep copy.
   */
  ThisType& operator=(const BackendType& other)
  {
    *backend_ = other.copy();
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
    return ThisType(backend_->copy());
  }

  void scal(const ScalarType& alpha)
  {
    const internal::VectorLockGuard DUNE_UNUSED(guard)(*mutexes_);
    backend() *= alpha;
  }

  void axpy(const ScalarType& alpha, const ThisType& xx)
  {
    if (xx.size() != size())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of x (" << xx.size() << ") does not match the size of this (" << size() << ")!");
    const internal::VectorLockGuard DUNE_UNUSED(guard)(*mutexes_);
    for (size_t ii = 0; ii < backend().size(); ++ii)
      backend()[ii] += alpha * xx.backend()[ii];
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

  inline void resize(const size_t new_size)
  {
    backend_->resize(new_size);
  }

  void add_to_entry(const size_t ii, const ScalarType& value)
  {
    assert(ii < size());
    internal::LockGuard DUNE_UNUSED(lock)(*mutexes_, ii, size());
    backend()[ii] += value;
  }

  void set_entry(const size_t ii, const ScalarType& value)
  {
    assert(ii < size());
    backend()[ii] = value;
  }

  ScalarType get_entry(const size_t ii) const
  {
    assert(ii < size());
    return backend()[ii];
  }

protected:
  inline ScalarType& get_unchecked_ref(const size_t ii)
  {
    return backend_->operator[](ii);
  }

  inline const ScalarType& get_unchecked_ref(const size_t ii) const
  {
    return backend_->operator[](ii);
  }

public:
  inline ScalarType& operator[](const size_t ii)
  {
    return get_unchecked_ref(ii);
  }

  inline const ScalarType& operator[](const size_t ii) const
  {
    return get_unchecked_ref(ii);
  }

  /// \}
  /// \name These methods override default implementations from VectorInterface.
  /// \{

  virtual ScalarType dot(const ThisType& other) const override final
  {
    if (other.size() != size())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    return backend() * other.backend();
  } // ... dot(...)

  virtual RealType l1_norm() const override final
  {
    return backend().l1_norm();
  }

  virtual RealType l2_norm() const override final
  {
    return backend().l2_norm();
  }

  virtual RealType sup_norm() const override final
  {
    return backend().sup_norm();
  }

  virtual void iadd(const ThisType& other) override final
  {
    if (other.size() != size())
      if (other.size() != size())
        DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                   "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    const internal::VectorLockGuard DUNE_UNUSED(guard)(*mutexes_);
    backend() += other.backend();
  } // ... iadd(...)

  virtual void isub(const ThisType& other) override final
  {
    if (other.size() != size())
      if (other.size() != size())
        DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                   "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    const internal::VectorLockGuard DUNE_UNUSED(guard)(*mutexes_);
    backend() -= other.backend();
  } // ... isub(...)

  // without these using declarations, the free operator+/* function in xt/common/vector.hh is chosen instead of the
  // member function
  using InterfaceType::operator+;
  using InterfaceType::operator-;
  using InterfaceType::operator*;

private:
  friend class VectorInterface<internal::CommonDenseVectorTraits<ScalarType>, ScalarType>;

  std::shared_ptr<BackendType> backend_;
  std::unique_ptr<MutexesType> mutexes_;
}; // class CommonDenseVector

} // namespace LA
namespace Common {


template <class T>
struct VectorAbstraction<LA::CommonDenseVector<T>>
  : public LA::internal::VectorAbstractionBase<LA::CommonDenseVector<T>>
{};


} // namespace Common
} // namespace XT
} // namespace Dune


// begin: this is what we need for the lib
#if DUNE_XT_WITH_PYTHON_BINDINGS


extern template class Dune::XT::LA::CommonDenseVector<double>;


#endif // DUNE_XT_WITH_PYTHON_BINDINGS
// end: this is what we need for the lib


#endif // DUNE_XT_LA_CONTAINER_COMMON_VECTOR_DENSE_HH
