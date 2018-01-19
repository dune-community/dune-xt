// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_XT_LA_CONTAINER_COMMON_VECTOR_DENSE_HH
#define DUNE_XT_LA_CONTAINER_COMMON_VECTOR_DENSE_HH

#include <cmath>
#include <initializer_list>
#include <memory>
#include <type_traits>
#include <vector>
#include <complex>
#include <mutex>

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

template <class ScalarImp>
class CommonDenseMatrix;

namespace internal {


/// Traits for CommonDenseVector
template <class ScalarImp>
class CommonDenseVectorTraits
{
public:
  typedef typename Dune::FieldTraits<ScalarImp>::field_type ScalarType;
  typedef typename Dune::FieldTraits<ScalarImp>::real_type RealType;
  typedef ScalarType DataType;
  typedef CommonDenseVector<ScalarType> derived_type;
  typedef Dune::DynamicVector<ScalarType> BackendType;
  static const Backends backend_type = Backends::common_dense;
  static const Backends dense_matrix_type = Backends::common_dense;
  static const Backends sparse_matrix_type = Backends::common_sparse;
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
  typedef typename Traits::DataType DataType;
  typedef typename Traits::BackendType BackendType;

  explicit CommonDenseVector(const size_t ss = 0, const ScalarType value = ScalarType(0), const size_t num_mutexes = 1)
    : backend_(new BackendType(ss, value))
    , mutexes_(num_mutexes > 0 ? std::make_shared<std::vector<std::mutex>>(num_mutexes) : nullptr)
    , unshareable_(false)
  {
  }

  explicit CommonDenseVector(const std::vector<ScalarType>& other, const size_t num_mutexes = 1)
    : backend_(new BackendType(other.size()))
    , mutexes_(num_mutexes > 0 ? std::make_shared<std::vector<std::mutex>>(num_mutexes) : nullptr)
    , unshareable_(false)
  {
    for (size_t ii = 0; ii < other.size(); ++ii)
      backend_->operator[](ii) = other[ii];
  }

  explicit CommonDenseVector(const std::initializer_list<ScalarType>& other, const size_t num_mutexes = 1)
    : backend_(new BackendType(other.size()))
    , mutexes_(num_mutexes > 0 ? std::make_shared<std::vector<std::mutex>>(num_mutexes) : nullptr)
    , unshareable_(false)
  {
    size_t ii = 0;
    for (auto element : other) {
      backend_->operator[](ii) = element;
      ++ii;
    }
  } // CommonDenseVector(...)

  CommonDenseVector(const ThisType& other)
    : backend_(other.unshareable_ ? std::make_shared<BackendType>(*other.backend_) : other.backend_)
    , mutexes_(other.unshareable_ ? std::make_shared<std::vector<std::mutex>>(other.mutexes_->size()) : other.mutexes_)
    , unshareable_(false)
  {
  }

  explicit CommonDenseVector(const BackendType& other,
                             const bool /*prune*/ = false,
                             const ScalarType /*eps*/ = Common::FloatCmp::DefaultEpsilon<ScalarType>::value(),
                             const size_t num_mutexes = 1)
    : backend_(new BackendType(other))
    , mutexes_(num_mutexes > 0 ? std::make_shared<std::vector<std::mutex>>(num_mutexes) : nullptr)
    , unshareable_(false)
  {
  }

  /**
   *  \note Takes ownership of backend_ptr in the sense that you must not delete it afterwards!
   */
  explicit CommonDenseVector(BackendType* backend_ptr, const size_t num_mutexes = 1)
    : backend_(backend_ptr)
    , mutexes_(num_mutexes > 0 ? std::make_shared<std::vector<std::mutex>>(num_mutexes) : nullptr)
    , unshareable_(false)
  {
  }

  explicit CommonDenseVector(std::shared_ptr<BackendType> backend_ptr, const size_t num_mutexes = 1)
    : backend_(backend_ptr)
    , mutexes_(num_mutexes > 0 ? std::make_shared<std::vector<std::mutex>>(num_mutexes) : nullptr)
    , unshareable_(false)
  {
  }

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      backend_ = other.unshareable_ ? std::make_shared<BackendType>(*other.backend_) : other.backend_;
      mutexes_ = other.mutexes_
                     ? (other.unshareable_ ? std::make_shared<std::vector<std::mutex>>(other.mutexes_->size())
                                           : other.mutexes_)
                     : nullptr;
      unshareable_ = false;
    }
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
    if (xx.size() != size())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of x (" << xx.size() << ") does not match the size of this (" << size() << ")!");
    for (size_t ii = 0; ii < backend_ref.size(); ++ii)
      backend_ref[ii] += alpha * xx.backend()[ii];
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
    auto& backend_ref = backend();
    internal::LockGuard DUNE_UNUSED(lock)(mutexes_, ii);
    assert(ii < size());
    backend_ref[ii] += value;
  }

  void set_entry(const size_t ii, const ScalarType& value)
  {
    auto& backend_ref = backend();
    assert(ii < size());
    backend_ref[ii] = value;
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
    ensure_uniqueness();
    unshareable_ = true;
    return get_unchecked_ref(ii);
  }

  inline const ScalarType& operator[](const size_t ii) const
  {
    ensure_uniqueness();
    unshareable_ = true;
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
    return backend().one_norm();
  }

  virtual RealType l2_norm() const override final
  {
    return backend().two_norm();
  }

  virtual RealType sup_norm() const override final
  {
    return backend().infinity_norm();
  }

  virtual void iadd(const ThisType& other) override final
  {
    auto& backend_ref = backend();
    const internal::VectorLockGuard DUNE_UNUSED(guard)(mutexes_);
    if (other.size() != size())
      if (other.size() != size())
        DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                   "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    backend_ref += other.backend();
  } // ... iadd(...)

  virtual void isub(const ThisType& other) override final
  {
    auto& backend_ref = backend();
    const internal::VectorLockGuard DUNE_UNUSED(guard)(mutexes_);
    if (other.size() != size())
      if (other.size() != size())
        DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                   "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    backend_ref -= other.backend();
  } // ... isub(...)

  // without these using declarations, the free operator+/* function in xt/common/vector.hh is chosen instead of the
  // member function
  using VectorInterfaceType::operator+;
  using VectorInterfaceType::operator-;
  using VectorInterfaceType::operator*;

protected:
  /**
   * \see ContainerInterface
   */
  inline void ensure_uniqueness() const
  {
    if (!backend_.unique()) {
      assert(!unshareable_);
      const internal::VectorLockGuard DUNE_UNUSED(guard)(mutexes_);
      if (!backend_.unique()) {
        backend_ = std::make_shared<BackendType>(*backend_);
        mutexes_ = mutexes_ ? std::make_shared<std::vector<std::mutex>>(mutexes_->size()) : nullptr;
      }
    }
  } // ... ensure_uniqueness(...)

private:
  friend class VectorInterface<internal::CommonDenseVectorTraits<ScalarType>, ScalarType>;
  friend class CommonDenseMatrix<ScalarType>;

  mutable std::shared_ptr<BackendType> backend_;
  mutable std::shared_ptr<std::vector<std::mutex>> mutexes_;
  mutable bool unshareable_;
}; // class CommonDenseVector

} // namespace LA
namespace Common {


template <class T>
struct VectorAbstraction<LA::CommonDenseVector<T>>
    : public LA::internal::VectorAbstractionBase<LA::CommonDenseVector<T>>
{
};


} // namespace Common
} // namespace XT
} // namespace Dune


// begin: this is what we need for the lib
#if DUNE_XT_WITH_PYTHON_BINDINGS


extern template class Dune::XT::LA::CommonDenseVector<double>;


#endif // DUNE_XT_WITH_PYTHON_BINDINGS
// end: this is what we need for the lib


#endif // DUNE_XT_LA_CONTAINER_COMMON_VECTOR_DENSE_HH
