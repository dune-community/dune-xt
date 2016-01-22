// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Barbara Verfürth (2015)
//   Felix Schindler  (2014 - 2015)
//   Rene Milk        (2015)
//   Tobias Leibner   (2014)

#ifndef DUNE_STUFF_LA_CONTAINER_EIGEN_BASE_HH
#define DUNE_STUFF_LA_CONTAINER_EIGEN_BASE_HH

#include <memory>
#include <type_traits>
#include <vector>
#include <complex>

#include <dune/stuff/common/disable_warnings.hh>
#if HAVE_EIGEN
#include <Eigen/Core>
#endif
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/common/typetraits.hh>
#include <dune/common/ftraits.hh>

#include <dune/stuff/aliases.hh>
#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/crtp.hh>
#include <dune/stuff/common/float_cmp.hh>

#include "dune/stuff/la/container/interfaces.hh"

namespace Dune {
namespace Stuff {
namespace LA {

template <class ScalarImp>
class EigenDenseMatrix;

template <class ScalarType>
class EigenRowMajorSparseMatrix;

#if HAVE_EIGEN

template <class Traits, class ScalarImp>
class EigenBaseVector;

/**
 *  \brief Base class for all eigen implementations of VectorInterface.
 */
template <class ImpTraits, class ScalarImp = double>
class EigenBaseVector : public VectorInterface<ImpTraits, ScalarImp>, public ProvidesBackend<ImpTraits>
{
  typedef VectorInterface<ImpTraits, ScalarImp> VectorInterfaceType;
  typedef EigenBaseVector<ImpTraits, ScalarImp> ThisType;

public:
  typedef ImpTraits Traits;
  typedef typename Traits::ScalarType ScalarType;
  typedef typename Traits::RealType RealType;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::derived_type VectorImpType;

  VectorImpType& operator=(const ThisType& other)
  {
    if (this != &other)
      backend_ = other.backend_;
    return this->as_imp();
  } // ... operator=(...)

  VectorImpType& operator=(const ScalarType& value)
  {
    ensure_uniqueness();
    for (auto& element : *this)
      element = value;
    return this->as_imp();
  } // ... operator=(...)

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

  VectorImpType copy() const
  {
    return VectorImpType(*backend_);
  }

  void scal(const ScalarType& alpha)
  {
    backend() *= alpha;
  }

  template <class T>
  void axpy(const ScalarType& alpha, const EigenBaseVector<T, ScalarType>& xx)
  {
    if (xx.size() != size())
      DUNE_THROW(Exceptions::shapes_do_not_match,
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
    return backend_->operator[](ii);
  }

protected:
  inline ScalarType& get_entry_ref(const size_t ii)
  {
    return backend()[ii];
  }

  inline const ScalarType& get_entry_ref(const size_t ii) const
  {
    return backend()[ii];
  }

  /// \}

public:
  /// \name These methods override default implementations from VectorInterface.
  /// \{

  virtual std::pair<size_t, RealType> amax() const override final
  {
    auto result            = std::make_pair(size_t(0), RealType(0));
    size_t min_index       = 0;
    size_t max_index       = 0;
    const RealType minimum = (backend_->cwiseAbs()).minCoeff(&min_index);
    const RealType maximum = (backend_->cwiseAbs()).maxCoeff(&max_index);
    if (maximum < minimum || (Common::FloatCmp::eq(maximum, minimum) && max_index > min_index)) {
      result.first  = min_index;
      result.second = minimum;
    } else {
      result.first  = max_index;
      result.second = maximum;
    }
    return result;
  } // ... amax(...)

  template <class T>
  ScalarType dot(const EigenBaseVector<T, ScalarType>& other) const
  {
    if (other.size() != size())
      DUNE_THROW(Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    return backend_->transpose() * *(other.backend_);
  } // ... dot(...)

  virtual ScalarType dot(const VectorImpType& other) const override final
  {
    return this->template dot<Traits>(other);
  }

  virtual RealType l1_norm() const override final
  {
    return backend_->template lpNorm<1>();
  }

  virtual RealType l2_norm() const override final
  {
    return backend_->template lpNorm<2>();
  }

  virtual RealType sup_norm() const override final
  {
    return backend_->template lpNorm<::Eigen::Infinity>();
  }

  template <class T1, class T2>
  void add(const EigenBaseVector<T1, ScalarType>& other, EigenBaseVector<T2, ScalarType>& result) const
  {
    if (other.size() != size())
      DUNE_THROW(Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    if (result.size() != size())
      DUNE_THROW(Exceptions::shapes_do_not_match,
                 "The size of result (" << result.size() << ") does not match the size of this (" << size() << ")!");
    result.backend() = *backend_ + *(other.backend_);
  } // ... add(...)

  virtual void add(const VectorImpType& other, VectorImpType& result) const override final
  {
    return this->template add<Traits, Traits>(other, result);
  }

  template <class T>
  void iadd(const EigenBaseVector<T, ScalarType>& other)
  {
    if (other.size() != size())
      DUNE_THROW(Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    backend() += *(other.backend_);
  } // ... iadd(...)

  virtual void iadd(const VectorImpType& other) override final
  {
    return this->template iadd<Traits>(other);
  }

  template <class T1, class T2>
  void sub(const EigenBaseVector<T1, ScalarType>& other, EigenBaseVector<T2, ScalarType>& result) const
  {
    if (other.size() != size())
      DUNE_THROW(Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    if (result.size() != size())
      DUNE_THROW(Exceptions::shapes_do_not_match,
                 "The size of result (" << result.size() << ") does not match the size of this (" << size() << ")!");
    result.backend() = *backend_ - *(other.backend_);
  } // ... sub(...)

  virtual void sub(const VectorImpType& other, VectorImpType& result) const override final
  {
    return this->template sub<Traits, Traits>(other, result);
  }

  template <class T>
  void isub(const EigenBaseVector<T, ScalarType>& other)
  {
    if (other.size() != size())
      DUNE_THROW(Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    backend() -= *(other.backend_);
  } // ... isub(...)

  virtual void isub(const VectorImpType& other) override final
  {
    this->template isub<Traits>(other);
  }

  /// \{

  //! disambiguation necessary since it exists in multiple bases
  using VectorInterfaceType::as_imp;

private:
  /**
   * \see ContainerInterface
   */
  void ensure_uniqueness() const
  {
    CHECK_AND_CALL_CRTP(VectorInterfaceType::as_imp().ensure_uniqueness());
    VectorInterfaceType::as_imp().ensure_uniqueness();
  }

#ifndef NDEBUG
  //! disambiguation necessary since it exists in multiple bases
  using VectorInterfaceType::crtp_mutex_;
#endif

  friend class VectorInterface<Traits, ScalarType>;
  friend class EigenDenseMatrix<ScalarType>;
  friend class EigenRowMajorSparseMatrix<ScalarType>;

protected:
  mutable std::shared_ptr<BackendType> backend_;
}; // class EigenBaseVector

#else // HAVE_EIGEN

template <class Traits, class ScalarImp>
class EigenBaseVector
{
  static_assert(Dune::AlwaysFalse<ScalarImp>::value, "You are missing Eigen!");
};

#endif // HAVE_EIGEN

} // namespace LA
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_LA_CONTAINER_EIGEN_BASE_HH
