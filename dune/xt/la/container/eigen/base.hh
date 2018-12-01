// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Barbara Verfürth (2015)
//   Felix Schindler  (2014 - 2017)
//   Rene Milk        (2015 - 2016, 2018)
//   Tobias Leibner   (2014, 2017 - 2018)

#ifndef DUNE_XT_LA_CONTAINER_EIGEN_BASE_HH
#define DUNE_XT_LA_CONTAINER_EIGEN_BASE_HH

#include <memory>
#include <type_traits>
#include <vector>
#include <complex>
#include <mutex>

#if HAVE_EIGEN
#include <dune/xt/common/disable_warnings.hh>
#include <Eigen/Core>
#include <dune/xt/common/reenable_warnings.hh>
#endif

#include <dune/common/typetraits.hh>
#include <dune/common/ftraits.hh>
#include <dune/common/unused.hh>

#include <dune/xt/common/crtp.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/float_cmp.hh>

#include "dune/xt/la/container/interfaces.hh"

namespace Dune {
namespace XT {
namespace LA {


template <class ScalarImp>
class EigenDenseMatrix;

template <class ScalarType>
class EigenRowMajorSparseMatrix;

#if HAVE_EIGEN


/**
 *  \brief Base class for all eigen implementations of VectorInterface.
 */
template <class ImpTraits, class ScalarImp = double>
class EigenBaseVector
  : public VectorInterface<ImpTraits, ScalarImp>
  , public ProvidesBackend<ImpTraits>
{

  using ThisType = EigenBaseVector;
  using InterfaceType = VectorInterface<ImpTraits, ScalarImp>;

public:
  using typename InterfaceType::RealType;
  using typename InterfaceType::ScalarType;
  using typename InterfaceType::Traits;
  using typename ProvidesBackend<Traits>::BackendType;
  using VectorImpType = typename InterfaceType::derived_type;

private:
  using MutexesType = typename Traits::MutexesType;

public:
  EigenBaseVector(size_t num_mutexes = 1)
    : mutexes_(std::make_unique<MutexesType>(num_mutexes))
  {}

  EigenBaseVector(const EigenBaseVector& other)
    : backend_(std::make_shared<BackendType>(*other.backend_))
    , mutexes_(std::make_unique<MutexesType>(other.mutexes_->size()))
  {}

  EigenBaseVector(EigenBaseVector&& source) = default;

  VectorImpType& operator=(const ThisType& other)
  {
    if (this != &other) {
      *backend_ = *other.backend_;
      mutexes_ = std::make_unique<MutexesType>(other.mutexes_->size());
    }
    return this->as_imp();
  } // ... operator=(...)

  VectorImpType& operator=(const ScalarType& value)
  {
    for (auto& element : *this)
      element = value;
    return this->as_imp();
  } // ... operator=(...)

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

  VectorImpType copy() const
  {
    return VectorImpType(*backend_);
  }

  void scal(const ScalarType& alpha)
  {
    const internal::VectorLockGuard DUNE_UNUSED(guard)(*mutexes_);
    backend() *= alpha;
  }

  template <class T>
  void axpy(const ScalarType& alpha, const EigenBaseVector<T, ScalarType>& xx)
  {
    const internal::VectorLockGuard DUNE_UNUSED(guard)(*mutexes_);
    if (xx.size() != size())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of xx (" << xx.size() << ") does not match the size of this (" << size() << ")!");
    backend() += alpha * xx.backend();
  } // ... axpy(...)

  bool has_equal_shape(const VectorImpType& other) const
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
    internal::LockGuard DUNE_UNUSED(lock)(*mutexes_, ii, size());
    assert(ii < size());
    backend()(ii) += value;
  }

  void set_entry(const size_t ii, const ScalarType& value)
  {
    assert(ii < size());
    backend()(ii) = value;
  }

  ScalarType get_entry(const size_t ii) const
  {
    assert(ii < size());
    return backend()[ii];
  }

protected:
  inline ScalarType& get_unchecked_ref(const size_t ii)
  {
    return (*backend_)[ii];
  }

  inline const ScalarType& get_unchecked_ref(const size_t ii) const
  {
    return (*backend_)[ii];
  }

public:
  inline ScalarType& operator[](const size_t ii)
  {
    return backend()[ii];
  }

  inline const ScalarType& operator[](const size_t ii) const
  {
    return backend()[ii];
  }

  /// \}

  /// \name These methods override default implementations from VectorInterface.
  /// \{

  virtual std::pair<size_t, RealType> amax() const override final
  {
    auto result = std::make_pair(size_t(0), RealType(0));
    size_t min_index = 0;
    size_t max_index = 0;
    const RealType minimum = (backend_->cwiseAbs()).minCoeff(&min_index);
    const RealType maximum = (backend_->cwiseAbs()).maxCoeff(&max_index);
    if (maximum < minimum || (Common::FloatCmp::eq(maximum, minimum) && max_index > min_index)) {
      result.first = min_index;
      result.second = minimum;
    } else {
      result.first = max_index;
      result.second = maximum;
    }
    return result;
  } // ... amax(...)

  template <class T>
  ScalarType dot(const EigenBaseVector<T, ScalarType>& other) const
  {
    if (other.size() != size())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    return backend().transpose() * other.backend();
  } // ... dot(...)

  virtual ScalarType dot(const VectorImpType& other) const override final
  {
    return this->template dot<Traits>(other);
  }

  virtual RealType l1_norm() const override final
  {
    return backend().template lpNorm<1>();
  }

  virtual RealType l2_norm() const override final
  {
    return backend().template lpNorm<2>();
  }

  virtual RealType sup_norm() const override final
  {
    return backend().template lpNorm<::Eigen::Infinity>();
  }

  template <class T>
  void iadd(const EigenBaseVector<T, ScalarType>& other)
  {
    const internal::VectorLockGuard DUNE_UNUSED(guard)(*mutexes_);
    if (other.size() != size())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    backend() += other.backend();
  } // ... iadd(...)

  virtual void iadd(const VectorImpType& other) override final
  {
    return this->template iadd<Traits>(other);
  }

  template <class T>
  void isub(const EigenBaseVector<T, ScalarType>& other)
  {
    const internal::VectorLockGuard DUNE_UNUSED(guard)(*mutexes_);
    if (other.size() != size())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    backend() -= other.backend();
  } // ... isub(...)

  virtual void isub(const VectorImpType& other) override final
  {
    this->template isub<Traits>(other);
  }

  /// \{

  //! disambiguation necessary since it exists in multiple bases
  using InterfaceType::as_imp;

private:
#ifndef NDEBUG
  //! disambiguation necessary since it exists in multiple bases
  using InterfaceType::crtp_mutex_;
#endif

  friend class VectorInterface<Traits, ScalarType>;
  friend class EigenDenseMatrix<ScalarType>;
  friend class EigenRowMajorSparseMatrix<ScalarType>;

protected:
  std::shared_ptr<BackendType> backend_;
  std::unique_ptr<MutexesType> mutexes_;
}; // class EigenBaseVector

#else // HAVE_EIGEN

template <class Traits, class ScalarImp>
class EigenBaseVector
{
  static_assert(Dune::AlwaysFalse<ScalarImp>::value, "You are missing Eigen!");
};

#endif // HAVE_EIGEN

} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_CONTAINER_EIGEN_BASE_HH
