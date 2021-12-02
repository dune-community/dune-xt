// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018 - 2020)
//   Tobias Leibner  (2017 - 2020)

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

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/math.hh>
#include <dune/xt/common/numeric.hh>

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

  // Constructs ss elements with value value.
  CommonDenseVectorBackend(const size_t ss = 0, const ScalarType& value = ScalarType())
    : size_(ss)
    , values_vector_(size_, value)
    , values_ptr_(values_vector_.data())
  {}

  // Uses the provided buffer with size ss. The pointer values_ptr should point to the first element of the buffer.
  // The user has to ensure that the buffer has at least ss elements, we have no way to check it here.
  // Does not take ownership of the buffer and also does not copy the buffer, i.e.,
  // the user has to ensure that the buffer lives at least as long as this class,
  // and the user is reponsible for cleaning up the memory.
  CommonDenseVectorBackend(const size_t ss, ScalarType* values_ptr)
    : size_(ss)
    , values_ptr_(values_ptr)
  {}

  // The copy constructor performs a deep copy of the values, even if other was constructed from a values_ptr
  // If you do not want an actual copy, use the constructor above.
  CommonDenseVectorBackend(const ThisType& other)
    : size_(other.size_)
    , values_vector_(size_)
    , values_ptr_(values_vector_.data())
  {
    std::copy_n(other.values_ptr(), size_, values_ptr());
  }

  CommonDenseVectorBackend(ThisType&& other) noexcept
    : size_(std::move(other.size_))
    , values_vector_(std::move(other.values_vector_))
    // If other was created from a ptr, we also want to use that pointer to avoid copying. The following assumes
    // values_vector_.data() points to the same location as other.values_vector_.data() did before the move. This is not
    // strictly guaranteed by the standard as of now but fulfilled in all major implementations (see
    // https://stackoverflow.com/a/25348988 and LWG open issue 2321).
    , values_ptr_(std::move(other.values_ptr_))
  {}

  ThisType copy() const
  {
    return ThisType(*this);
  }

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      size_ = other.size_;
      values_vector_ = std::vector<ScalarType>(size_);
      std::copy_n(other.values_ptr(), size_, values_vector_.begin());
      values_ptr_ = values_vector_.data();
    }
    return *this;
  }

  ScalarType& operator[](const size_t ii)
  {
    return values_ptr()[ii];
  }

  const ScalarType& operator[](const size_t ii) const
  {
    return values_ptr()[ii];
  }

  void resize(const size_t ss)
  {
    DUNE_THROW_IF(values_ptr_ != values_vector_.data(),
                  Dune::InvalidStateException,
                  "Cannot resize vector that has been initialized from a pointer!");
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
      values_ptr()[ii] += other[ii];
    return *this;
  }

  ThisType& operator-=(const ThisType& other)
  {
    assert(other.size_ == size_);
    for (size_t ii = 0; ii < size_; ++ii)
      values_ptr()[ii] -= other[ii];
    return *this;
  }

  ThisType& operator*=(const ThisType& other)
  {
    assert(other.size_ == size_);
    for (size_t ii = 0; ii < size_; ++ii)
      values_ptr()[ii] *= other[ii];
    return *this;
  }

  ThisType& operator*=(const ScalarType& val)
  {
    for (size_t ii = 0; ii < size_; ++ii)
      values_ptr()[ii] *= val;
    return *this;
  }

  ThisType& operator/=(const ThisType& other)
  {
    assert(other.size_ == size_);
    for (size_t ii = 0; ii < size_; ++ii)
      values_ptr()[ii] /= other[ii];
    return *this;
  }

  ScalarType operator*(const ThisType& other) const
  {
    assert(other.size() == size_);
    return XT::Common::transform_reduce(values_ptr(), values_ptr() + size_, other.values_ptr(), ScalarType(0.));
  }

  RealType l1_norm() const
  {
    using Dune::XT::Common::abs;
    if constexpr (std::is_convertible_v<ScalarType, RealType>) {
      return Common::reduce(values_ptr(), values_ptr() + size_, RealType(0.), [](const auto& a, const auto& b) {
        return abs(a) + abs(b);
      });
    } else {
      return std::accumulate(values_ptr(),
                             values_ptr() + size_,
                             RealType(0.),
                             [](const RealType& a, const ScalarType& b) { return a + abs(b); });
    }
  }

  RealType l2_norm() const
  {
    using Dune::XT::Common::abs;
    return std::sqrt(
        std::accumulate(values_ptr(), values_ptr() + size_, RealType(0.), [](const RealType& a, const ScalarType& b) {
          return a + std::pow(abs(b), 2);
        }));
  }

  RealType sup_norm() const
  {
    using Dune::XT::Common::abs;
    return abs(*std::max_element(
        values_ptr(), values_ptr() + size_, [](const ScalarType& a, const ScalarType& b) { return abs(a) < abs(b); }));
  }

  inline ScalarType* values_ptr()
  {
    assert(values_ptr_);
    return values_ptr_;
  }

  inline const ScalarType* values_ptr() const
  {
    assert(values_ptr_);
    return values_ptr_;
  }

  size_t size_;
  std::vector<ScalarType> values_vector_;
  ScalarType* values_ptr_;
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

  // the bool is only to disambiguate the constructors
  explicit CommonDenseVector(const size_t ss, ScalarType* values_ptr, const size_t num_mutexes, bool /*dummy*/)
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
    assert(backend_);
    return *backend_;
  }

  const BackendType& backend() const
  {
    assert(backend_);
    return *backend_;
  }

  /// \}
  /// \name Required by ProvidesDataAccess.
  /// \{

  /** \attention This makes only sense for scalar data types, not for complex! **/
  ScalarType* data()
  {
    return backend().values_ptr();
  }

  const ScalarType* data() const
  {
    return backend().values_ptr();
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
    [[maybe_unused]] const internal::VectorLockGuard guard(*mutexes_);
    backend() *= alpha;
  }

  void axpy(const ScalarType& alpha, const ThisType& xx)
  {
    const size_t sz = size();
    if (xx.size() != sz)
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of x (" << xx.size() << ") does not match the size of this (" << size() << ")!");
    [[maybe_unused]] const internal::VectorLockGuard guard(*mutexes_);
    auto* data_ptr = data();
    const auto* other_data_ptr = xx.data();
    for (size_t ii = 0; ii < sz; ++ii)
      data_ptr[ii] += alpha * other_data_ptr[ii];
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
    [[maybe_unused]] const internal::LockGuard lock(*mutexes_, ii, size());
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

  ScalarType dot(const ThisType& other) const final
  {
    if (other.size() != size())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    return backend() * other.backend();
  } // ... dot(...)

  RealType l1_norm() const final
  {
    return backend().l1_norm();
  }

  RealType l2_norm() const final
  {
    return backend().l2_norm();
  }

  RealType sup_norm() const final
  {
    return backend().sup_norm();
  }

  void iadd(const ThisType& other) final
  {
    if (other.size() != size())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    [[maybe_unused]] const internal::VectorLockGuard guard(*mutexes_);
    backend() += other.backend();
  } // ... iadd(...)

  void isub(const ThisType& other) final
  {
    if (other.size() != size())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    [[maybe_unused]] const internal::VectorLockGuard guard(*mutexes_);
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
