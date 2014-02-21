// This file is part of the dune-stuff project:
//   https://users.dune-project.org/projects/dune-stuff/
// Copyright Holders: Felix Albrecht, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_LA_CONTAINER_EIGEN_HH
#define DUNE_STUFF_LA_CONTAINER_EIGEN_HH

#include <memory>
#include <type_traits>

#if HAVE_EIGEN
#include <Eigen/Core>
#include <Eigen/SparseCore>
#endif // HAVE_EIGEN

#include <dune/common/typetraits.hh>

#include <dune/stuff/aliases.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/crtp.hh>

#include "interfaces.hh"
#include "pattern.hh"

namespace Dune {
namespace Pymor {
namespace Operators {

// forwards, needed for friendlyness
template <class ScalarImp>
class EigenRowMajorSparseInverse;

template <class ScalarImp>
class EigenRowMajorSparse;

} // namespace Operators
} // namespace Pymor
namespace Stuff {
namespace LA {


// forwards
template <class VectorImpTraits, class ScalarImp>
class EigenBaseVector;

template <class ScalarImp>
class EigenDenseVector;

template <class T>
class EigenMappedDenseVector;

template <class ScalarImp>
class EigenDenseMatrix;

template <class ScalarType>
class EigenRowMajorSparseMatrix;


class EigenVectorInterfaceDynamic
{
};
class EigenMatrixInterfaceDynamic
{
};


#if HAVE_EIGEN


/**
 *  \brief Traits for EigenDenseVector.
 */
template <class VectorImpTraits, class ScalarImp = double>
class EigenBaseVectorTraits
{
  static_assert(std::is_same<ScalarImp, typename VectorImpTraits::ScalarType>::value, "Types do not match!");

public:
  typedef ScalarImp ScalarType;
  typedef EigenBaseVector<VectorImpTraits, ScalarType> derived_type;
  typedef typename VectorImpTraits::BackendType BackendType;
  typedef typename VectorImpTraits::derived_type VectorImpType;
}; // class EigenBaseVectorTraits


/**
 *  \brief A dense vector implementation of VectorInterface using the eigen backend.
 */
template <class VectorImpTraits, class ScalarImp = double>
class EigenBaseVector : public VectorInterface<EigenBaseVectorTraits<VectorImpTraits, ScalarImp>>,
                        public EigenVectorInterfaceDynamic,
                        public ProvidesBackend<EigenBaseVectorTraits<VectorImpTraits, ScalarImp>>,
                        protected CRTPInterface<EigenBaseVector<VectorImpTraits, ScalarImp>, VectorImpTraits>
{
  typedef EigenBaseVector<VectorImpTraits, ScalarImp> ThisType;
  typedef VectorInterface<EigenBaseVectorTraits<VectorImpTraits, ScalarImp>> VectorInterfaceType;

public:
  typedef EigenBaseVectorTraits<VectorImpTraits> Traits;
  typedef typename Traits::ScalarType ScalarType;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::VectorImpType VectorImpType;

protected:
  typedef CRTPInterface<EigenBaseVector<VectorImpTraits, ScalarImp>, VectorImpTraits> CRTP;

public:
  /**
   * \defgroup backend ´´These methods are required by the ProvidesBackend interface.``
   * \{
   */
  BackendType& backend()
  {
    CHECK_CRTP(CRTP::as_imp(*this).backend());
    return CRTP::as_imp(*this).backend();
  } // ... backend(...)

  const BackendType& backend() const
  {
    CHECK_CRTP(CRTP::as_imp(*this).backend());
    return CRTP::as_imp(*this).backend();
  } // ... backend(...)
  /**
   * \}
   */

  /**
   * \defgroup container ´´These methods are required by ContainerInterface.``
   * \{
   */
  VectorImpType copy() const
  {
    return VectorImpType(backend());
  }

  void scal(const ScalarType& alpha)
  {
    backend() *= alpha;
  }

  template <class T>
  void axpy(const ScalarType& alpha, const EigenBaseVector<T, ScalarType>& xx)
  {
    if (xx.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of xx (" << xx.size() << ") does not match the size of this (" << size() << ")!");
    backend() += alpha * xx.backend();
  } // ... axpy(...)

  bool has_equal_shape(const ThisType& other) const
  {
    return size() == other.size();
  }
  /**
   * \}
   */

  /**
   * \defgroup vector_required ´´These methods are required by VectorInterface.``
   * \{
   */
  inline size_t size() const
  {
    return backend().size();
  }

  void add_to_entry(const size_t ii, const ScalarType& value)
  {
    assert(ii < size());
    backend()(ii) += value;
  } // ... add_to_entry(...)

  void set_entry(const size_t ii, const ScalarType& value)
  {
    assert(ii < size());
    backend()(ii) = value;
  } // ... set_entry(...)

  ScalarType get_entry(const size_t ii) const
  {
    assert(ii < size());
    return backend()(ii);
  } // ... get_entry(...)

protected:
  inline ScalarType& get_entry_ref(const size_t ii)
  {
    return backend()(ii);
  }

  inline const ScalarType& get_entry_ref(const size_t ii) const
  {
    return backend()(ii);
  }
  /**
   * \}
   */

public:
  /**
   * \defgroup vector_overrides ´´These methods override default implementations from VectorInterface.``
   * \{
   */
  virtual std::pair<size_t, ScalarType> amax() const DS_OVERRIDE
  {
    auto result              = std::make_pair(size_t(0), ScalarType(0));
    size_t min_index         = 0;
    size_t max_index         = 0;
    const ScalarType minimum = backend().minCoeff(&min_index);
    const ScalarType maximum = backend().maxCoeff(&max_index);
    if (std::abs(maximum) >= std::abs(minimum)) {
      result.first  = min_index;
      result.second = std::abs(minimum);
    } else {
      result.first  = max_index;
      result.second = maximum;
    }
    return result;
  } // ... amax(...)

  template <class T>
  bool almost_equal(const EigenBaseVector<T, ScalarType>& other,
                    const ScalarType epsilon = Dune::FloatCmp::DefaultEpsilon<ScalarType>::value()) const
  {
    if (other.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of other (" << other.size() << ") does not match the size of this (" << size()
                                                  << ")!");
    for (size_t ii = 0; ii < size(); ++ii)
      if (!Dune::FloatCmp::eq<ScalarType>(get_entry(ii), other.get_entry(ii), epsilon))
        return false;
    return true;
  } // ... almost_equal(...)

  virtual bool
  almost_equal(const ThisType& other,
               const ScalarType epsilon = Dune::FloatCmp::DefaultEpsilon<ScalarType>::value()) const DS_OVERRIDE
  {
    return this->template almost_equal<VectorImpTraits>(other, epsilon);
  }

  using VectorInterfaceType::almost_equal;

  template <class T>
  ScalarType dot(const EigenBaseVector<T, ScalarType>& other) const
  {
    if (other.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of other (" << other.size() << ") does not match the size of this (" << size()
                                                  << ")!");
    return backend().transpose() * other.backend();
  } // ... dot(...)

  virtual ScalarType dot(const ThisType& other) const DS_OVERRIDE
  {
    return this->template dot<VectorImpTraits>(other);
  }

  virtual ScalarType l1_norm() const DS_OVERRIDE
  {
    return backend().template lpNorm<1>();
  }

  virtual ScalarType l2_norm() const DS_OVERRIDE
  {
    return backend().template lpNorm<2>();
  }

  virtual ScalarType sup_norm() const DS_OVERRIDE
  {
    return backend().template lpNorm<::Eigen::Infinity>();
  }

  template <class T1, class T2>
  void add(const EigenBaseVector<T1, ScalarType>& other, EigenBaseVector<T2, ScalarType>& result) const
  {
    if (other.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of other (" << other.size() << ") does not match the size of this (" << size()
                                                  << ")!");
    if (result.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of result (" << result.size() << ") does not match the size of this (" << size()
                                                   << ")!");
    result.backend() = backend() + other.backend();
  } // ... add(...)

  virtual void add(const ThisType& other, ThisType& result) const DS_OVERRIDE
  {
    return this->template add<VectorImpTraits, VectorImpTraits>(other, result);
  }

  template <class T>
  VectorImpType add(const EigenBaseVector<T, ScalarType>& other) const
  {
    if (other.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of other (" << other.size() << ") does not match the size of this (" << size()
                                                  << ")!");
    return VectorImpType(backend() + other.backend());
  } // ... add(...)

  virtual ThisType add(const ThisType& other) const DS_OVERRIDE
  {
    return this->template add<VectorImpTraits>(other);
  }

  template <class T>
  void iadd(const EigenBaseVector<T, ScalarType>& other)
  {
    if (other.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of other (" << other.size() << ") does not match the size of this (" << size()
                                                  << ")!");
    backend() += other.backend();
  } // ... iadd(...)

  virtual void iadd(const ThisType& other) DS_OVERRIDE
  {
    return this->template iadd<VectorImpTraits>(other);
  }

  template <class T1, class T2>
  void sub(const EigenBaseVector<T1, ScalarType>& other, EigenBaseVector<T2, ScalarType>& result) const
  {
    if (other.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of other (" << other.size() << ") does not match the size of this (" << size()
                                                  << ")!");
    if (result.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of result (" << result.size() << ") does not match the size of this (" << size()
                                                   << ")!");
    result.backend() = backend() - other.backend();
  } // ... sub(...)

  virtual void sub(const ThisType& other, ThisType& result) const DS_OVERRIDE
  {
    return this->template sub<VectorImpTraits, VectorImpTraits>(other, result);
  }

  template <class T>
  VectorImpType sub(const EigenBaseVector<T, ScalarType>& other) const
  {
    if (other.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of other (" << other.size() << ") does not match the size of this (" << size()
                                                  << ")!");
    return VectorImpType(backend() - other.backend());
  } // ... sub(...)

  virtual ThisType sub(const ThisType& other) const DS_OVERRIDE
  {
    return this->template sub<VectorImpTraits>(other);
  }

  template <class T>
  void isub(const EigenBaseVector<T, ScalarType>& other)
  {
    if (other.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of other (" << other.size() << ") does not match the size of this (" << size()
                                                  << ")!");
    backend() -= other.backend();
  } // ... isub(...)

  virtual void isub(const ThisType& other) DS_OVERRIDE
  {
    this->template isub<VectorImpTraits>(other);
  }
  /**
   * \}
   */

private:
  friend class VectorInterface<EigenBaseVectorTraits<VectorImpTraits, ScalarType>>;
}; // class EigenBaseVector


/**
 *  \brief Traits for EigenDenseVector.
 */
template <class ScalarImp = double>
class EigenDenseVectorTraits
{
public:
  typedef ScalarImp ScalarType;
  typedef EigenDenseVector<ScalarType> derived_type;
  typedef typename ::Eigen::Matrix<ScalarType, ::Eigen::Dynamic, 1> BackendType;
}; // class EigenDenseVectorTraits


/**
 *  \brief A dense vector implementation of VectorInterface using the eigen backend.
 */
template <class ScalarImp = double>
class EigenDenseVector : public VectorInterface<EigenDenseVectorTraits<ScalarImp>>,
                         public EigenBaseVector<EigenDenseVectorTraits<ScalarImp>>,
                         public ProvidesDataAccess<EigenDenseVectorTraits<ScalarImp>>
{
  typedef EigenDenseVector<ScalarImp> ThisType;
  typedef VectorInterface<EigenDenseVectorTraits<ScalarImp>> VectorInterfaceType;
  typedef EigenBaseVector<EigenDenseVectorTraits<ScalarImp>> BaseType;
  static_assert(!std::is_same<DUNE_STUFF_SSIZE_T, int>::value,
                "You have to manually disable the constructor below which uses DUNE_STUFF_SSIZE_T!");

public:
  typedef EigenDenseVectorTraits<ScalarImp> Traits;
  typedef typename Traits::ScalarType ScalarType;
  typedef typename Traits::BackendType BackendType;

  EigenDenseVector(const size_t ss = 0, const ScalarType value = ScalarType(0))
    : backend_(new BackendType(ss))
  {
    if (FloatCmp::eq(value, ScalarType(0)))
      backend_->setZero();
    else {
      backend_->setOnes();
      backend_->operator*=(value);
    }
  }

  /// This constructor is needed for the python bindings.
  EigenDenseVector(const DUNE_STUFF_SSIZE_T ss, const ScalarType value = ScalarType(0))
    : backend_(new BackendType(VectorInterfaceType::assert_is_size_t_compatible_and_convert(ss)))
  {
    if (FloatCmp::eq(value, ScalarType(0)))
      backend_->setZero();
    else {
      backend_->setOnes();
      backend_->operator*=(value);
    }
  }

  EigenDenseVector(const int ss, const ScalarType value = ScalarType(0))
    : backend_(new BackendType(VectorInterfaceType::assert_is_size_t_compatible_and_convert(ss)))
  {
    if (FloatCmp::eq(value, ScalarType(0)))
      backend_->setZero();
    else {
      backend_->setOnes();
      backend_->operator*=(value);
    }
  }

  EigenDenseVector(const ThisType& other)
    : backend_(other.backend_)
  {
  }

  EigenDenseVector(const BackendType& other)
    : backend_(new BackendType(other))
  {
  }

  /**
   *  \note Takes ownership of backend_ptr in the sense that you must not delete it afterwards!
   */
  EigenDenseVector(BackendType* backend_ptr)
    : backend_(backend_ptr)
  {
  }

  EigenDenseVector(std::shared_ptr<BackendType> backend_ptr)
    : backend_(backend_ptr)
  {
  }

  ThisType& operator=(const ThisType& other)
  {
    backend_ = other.backend_;
    return *this;
  } // ... operator=(...)

  /**
   *  \note Does a deep copy.
   */
  ThisType& operator=(const BackendType& other)
  {
    backend_ = std::make_shared<BackendType>(other);
    return *this;
  } // ... operator=(...)

  /**
   * \defgroup backend ´´These methods are required by the ProvidesBackend interface.``
   * \{
   */
  BackendType& backend()
  {
    ensure_uniqueness();
    return *backend_;
  } // ... backend(...)

  const BackendType& backend() const
  {
    const_cast<ThisType&>(*this).ensure_uniqueness();
    return *backend_;
  } // ... backend(...)
  /**
   * \}
   */

  /**
   * \defgroup data ´´These methods are required by the ProvidesDataAccess interface.``
   * \{
   */
  ScalarType* data()
  {
    return backend_->data();
  }
  /**
   * \}
   */

  /**
   * \defgroup vector_required ´´These methods are required by VectorInterface.``
   * \{
   */
  using BaseType::copy;
  using BaseType::scal;
  /**
   * \}
   */

  /**
   * \defgroup vector_required ´´These methods are required by VectorInterface.``
   * \{
   */
  using BaseType::size;
  using BaseType::add_to_entry;
  using BaseType::set_entry;
  using BaseType::get_entry;

protected:
  using BaseType::get_entry_ref;
  /**
   * \}
   */

public:
  /**
   * \defgroup vector_overrides ´´These methods override default implementations from VectorInterface.``
   * \{
   */
  virtual bool
  almost_equal(const ThisType& other,
               const ScalarType epsilon = Dune::FloatCmp::DefaultEpsilon<ScalarType>::value()) const DS_OVERRIDE
  {
    return BaseType::template almost_equal<Traits>(other, epsilon);
  }

  using BaseType::almost_equal;

  virtual ScalarType dot(const ThisType& other) const DS_OVERRIDE
  {
    BaseType::template dot<Traits>(other);
  }

  using BaseType::dot;

  virtual void add(const ThisType& other, ThisType& result) const DS_OVERRIDE
  {
    BaseType::template add<Traits, Traits>(other, result);
  }

  using BaseType::sup_norm;

  virtual ThisType add(const ThisType& other) const DS_OVERRIDE
  {
    return BaseType::template add<Traits>(other);
  }

  using BaseType::add;

  virtual void iadd(const ThisType& other) DS_OVERRIDE
  {
    BaseType::template iadd<Traits>(other);
  }

  using BaseType::iadd;

  virtual void sub(const ThisType& other, ThisType& result) const DS_OVERRIDE
  {
    BaseType::template sub<Traits, Traits>(other, result);
  }

  virtual ThisType sub(const ThisType& other) const DS_OVERRIDE
  {
    return BaseType::template sub<Traits>(other);
  }

  using BaseType::sub;

  virtual void isub(const ThisType& other) DS_OVERRIDE
  {
    BaseType::template isub<Traits>(other);
  }

  using BaseType::isub;
  /**
   * \}
   */

  using BaseType::operator*=;

private:
  inline void ensure_uniqueness()
  {
    if (!backend_.unique())
      backend_ = std::make_shared<BackendType>(*backend_);
  } // ... ensure_uniqueness(...)

  friend class VectorInterface<EigenDenseVectorTraits<ScalarType>>;
  friend class EigenBaseVector<EigenDenseVectorTraits<ScalarType>>;

  std::shared_ptr<BackendType> backend_;
}; // class EigenDenseVector


/**
 *  \brief Traits for EigenMappedDenseVector.
 */
template <class ScalarImp = double>
class EigenMappedDenseVectorTraits
{
  typedef typename ::Eigen::Matrix<ScalarImp, ::Eigen::Dynamic, 1> PlainBackendType;

public:
  typedef ScalarImp ScalarType;
  typedef EigenMappedDenseVector<ScalarType> derived_type;
  typedef Eigen::Map<PlainBackendType> BackendType;
};


/**
 *  \brief  A dense vector implementation of VectorInterface using the eigen backend which wrappes a raw array.
 */
template <class ScalarImp = double>
class EigenMappedDenseVector : public VectorInterface<EigenMappedDenseVectorTraits<ScalarImp>>,
                               public EigenVectorInterfaceDynamic,
                               public ProvidesBackend<EigenMappedDenseVectorTraits<ScalarImp>>
{
  typedef EigenMappedDenseVector<ScalarImp> ThisType;
  typedef VectorInterface<EigenMappedDenseVectorTraits<ScalarImp>> VectorInterfaceType;
  static_assert(std::is_same<ScalarImp, double>::value, "Undefined behaviour for non-double data!");
  static_assert(!std::is_same<DUNE_STUFF_SSIZE_T, int>::value,
                "You have to manually disable the constructor below which uses DUNE_STUFF_SSIZE_T!");

public:
  typedef EigenMappedDenseVectorTraits<ScalarImp> Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::ScalarType ScalarType;

  /**
   *  \brief  This is the constructor of interest which wrappes a raw array.
   */
  EigenMappedDenseVector(ScalarType* data, size_t data_size)
    : backend_(new BackendType(data, data_size))
  {
  }

  /**
   *  \brief  This constructor allows to create an instance of this type just like any other vector.
   */
  EigenMappedDenseVector(const size_t ss = 0, const ScalarType value = ScalarType(0))
    : backend_(new BackendType(new ScalarType[ss], ss))
  {
    if (FloatCmp::eq(value, ScalarType(0)))
      backend_->setZero();
    else {
      backend_->setOnes();
      backend_->operator*=(value);
    }
  }

  /// This constructor is needed for the python bindings.
  EigenMappedDenseVector(const DUNE_STUFF_SSIZE_T ss, const ScalarType value = ScalarType(0))
    : backend_(new BackendType(new ScalarType[VectorInterfaceType::assert_is_size_t_compatible_and_convert(ss)],
                               VectorInterfaceType::assert_is_size_t_compatible_and_convert(ss)))
  {
    if (FloatCmp::eq(value, ScalarType(0)))
      backend_->setZero();
    else {
      backend_->setOnes();
      backend_->operator*=(value);
    }
  }

  EigenMappedDenseVector(const int ss, const ScalarType value = ScalarType(0))
    : backend_(new BackendType(new ScalarType[VectorInterfaceType::assert_is_size_t_compatible_and_convert(ss)],
                               VectorInterfaceType::assert_is_size_t_compatible_and_convert(ss)))
  {
    if (FloatCmp::eq(value, ScalarType(0)))
      backend_->setZero();
    else {
      backend_->setOnes();
      backend_->operator*=(value);
    }
  }


  /**
   *  \brief  This constructor does not do a deep copy.
   */
  EigenMappedDenseVector(const ThisType& other)
    : backend_(other.backend_)
  {
  }

  /**
   * \brief This constructor does a deep copy.
   */
  EigenMappedDenseVector(const BackendType& other)
    : backend_(new BackendType(new ScalarType[other.size()], other.size()))
  {
    backend_->operator=(other);
  }

  /**
   *  \note Takes ownership of backend_ptr in the sense that you must not delete it afterwards!
   */
  EigenMappedDenseVector(BackendType* backend_ptr)
    : backend_(backend_ptr)
  {
  }

  EigenMappedDenseVector(std::shared_ptr<BackendType> backend_ptr)
    : backend_(backend_ptr)
  {
  }

  /**
   *  \brief  This does not do a deep copy.
   */
  ThisType& operator=(const ThisType& other)
  {
    backend_ = other.backend_;
    return *this;
  } // ... operator=(...)

  /**
   * \defgroup backend ´´These methods are required by the ProvidesBackend interface.``
   * \{
   */

  BackendType& backend()
  {
    ensure_uniqueness();
    return *backend_;
  } // ... backend(...)

  const BackendType& backend() const
  {
    const_cast<ThisType&>(*this).ensure_uniqueness();
    return *backend_;
  } // ... backend(...)

  /**
   * \}
   */

  /**
   * \defgroup container ´´These methods are required by ContainerInterface.``
   * \{
   */
  ThisType copy() const
  {
    return ThisType(*backend_);
  }

  void scal(const ScalarType& alpha)
  {
    ensure_uniqueness();
    backend_->operator*=(alpha);
  } // ... scal(...)

  void axpy(const ScalarType& alpha, const ThisType& xx)
  {
    if (xx.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of xx (" << xx.size() << ") does not match the size of this (" << size() << ")!");
    ensure_uniqueness();
    auto& this_ref     = *backend_;
    const auto& xx_ref = *(xx.backend_);
    this_ref += alpha * xx_ref;
  } // ... axpy(...)

  bool has_equal_shape(const ThisType& other) const
  {
    return size() == other.size();
  }
  /**
   * \}
   */

  /**
   * \defgroup vector_required ´´These methods are required by VectorInterface.``
   * \{
   */

  inline size_t size() const
  {
    return backend_->size();
  }

  void add_to_entry(const size_t ii, const ScalarType& value)
  {
    assert(ii < size());
    ensure_uniqueness();
    backend_->operator()(ii) += value;
  } // ... add_to_entry(...)

  void set_entry(const size_t ii, const ScalarType& value)
  {
    assert(ii < size());
    ensure_uniqueness();
    backend_->operator()(ii) = value;
  } // ... set_entry(...)

  ScalarType get_entry(const size_t ii) const
  {
    assert(ii < size());
    return backend_->operator()(ii);
  } // ... get_entry(...)

private:
  inline ScalarType& get_entry_ref(const size_t ii)
  {
    return backend_->operator()(ii);
  }

  inline const ScalarType& get_entry_ref(const size_t ii) const
  {
    return backend_->operator()(ii);
  }

public:
  /**
   * \}
   */

  /**
   * \defgroup vector_overrides ´´These methods override default implementations from VectorInterface.``
   * \{
   */

  virtual std::pair<size_t, ScalarType> amax() const DS_OVERRIDE
  {
    auto result              = std::make_pair(size_t(0), ScalarType(0));
    size_t min_index         = 0;
    size_t max_index         = 0;
    const ScalarType minimum = backend_->minCoeff(&min_index);
    const ScalarType maximum = backend_->maxCoeff(&max_index);
    if (std::abs(maximum) >= std::abs(minimum)) {
      result.first  = min_index;
      result.second = std::abs(minimum);
    } else {
      result.first  = max_index;
      result.second = maximum;
    }
    return result;
  } // ... amax(...)

  bool almost_equal(const EigenDenseVector<ScalarType>& other,
                    const ScalarType epsilon = Dune::FloatCmp::DefaultEpsilon<ScalarType>::value()) const
  {
    if (other.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of other (" << other.size() << ") does not match the size of this (" << size()
                                                  << ")!");
    for (size_t ii = 0; ii < size(); ++ii)
      if (!Dune::FloatCmp::eq<ScalarType>(backend_->operator()(ii), other.backend_->operator()(ii), epsilon))
        return false;
    return true;
  } // ... almost_equal(...)

  using VectorInterfaceType::almost_equal;

  virtual ScalarType dot(const ThisType& other) const DS_OVERRIDE
  {
    if (other.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of other (" << other.size() << ") does not match the size of this (" << size()
                                                  << ")!");
    const auto& other_ref = *(other.backend_);
    return backend_->transpose() * other_ref;
  } // ... dot(...)

  virtual ScalarType l1_norm() const DS_OVERRIDE
  {
    return backend_->template lpNorm<1>();
  }

  virtual ScalarType l2_norm() const DS_OVERRIDE
  {
    return backend_->template lpNorm<2>();
  }

  virtual ScalarType sup_norm() const DS_OVERRIDE
  {
    return backend_->template lpNorm<::Eigen::Infinity>();
  }

  virtual void add(const ThisType& other, ThisType& result) const DS_OVERRIDE
  {
    if (other.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of other (" << other.size() << ") does not match the size of this (" << size()
                                                  << ")!");
    if (result.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of result (" << result.size() << ") does not match the size of this (" << size()
                                                   << ")!");
    result.backend() = backend_->operator+(*(other.backend_));
  } // ... add(...)

  virtual void iadd(const ThisType& other) DS_OVERRIDE
  {
    if (other.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of other (" << other.size() << ") does not match the size of this (" << size()
                                                  << ")!");
    ensure_uniqueness();
    const auto& other_ref = *(other.backend_);
    backend_->operator+=(other_ref);
  } // ... iadd(...)

  virtual void sub(const ThisType& other, ThisType& result) const DS_OVERRIDE
  {
    if (other.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of other (" << other.size() << ") does not match the size of this (" << size()
                                                  << ")!");
    if (result.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of result (" << result.size() << ") does not match the size of this (" << size()
                                                   << ")!");
    result.backend() = backend_->operator-(*(other.backend_));
  } // ... sub(...)

  void isub(const ThisType& other) DS_OVERRIDE
  {
    if (other.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of other (" << other.size() << ") does not match the size of this (" << size()
                                                  << ")!");
    ensure_uniqueness();
    const auto& other_ref = *(other.backend_);
    backend_->operator-=(other_ref);
  } // ... isub(...)

  /**
   * \}
   */
private:
  inline void ensure_uniqueness()
  {
    if (!backend_.unique()) {
      auto new_backend     = std::make_shared<BackendType>(new ScalarType[size()], size());
      new_backend->operator=(*backend_);
      backend_             = new_backend;
    }
  } // ... ensure_uniqueness(...)

  friend class VectorInterface<EigenMappedDenseVectorTraits<ScalarType>>;
  friend class EigenDenseVector<ScalarType>;
  friend class EigenDenseMatrix<ScalarType>;
  friend class EigenRowMajorSparseMatrix<ScalarType>;
  friend class Dune::Pymor::Operators::EigenRowMajorSparseInverse<ScalarType>;
  friend class Dune::Pymor::Operators::EigenRowMajorSparse<ScalarType>;

  std::shared_ptr<BackendType> backend_;
}; // class EigenMappedDenseVector


/**
 *  \brief Traits for EigenDenseMatrix.
 */
template <class ScalarImp = double>
class EigenDenseMatrixTraits
{
public:
  typedef ScalarImp ScalarType;
  typedef EigenDenseMatrix<ScalarType> derived_type;
  typedef typename ::Eigen::Matrix<ScalarType, ::Eigen::Dynamic, ::Eigen::Dynamic> BackendType;
}; // class DenseMatrixTraits


/**
 *  \brief  A dense metrix implementation of MatrixInterface using the eigen backend.
 */
template <class ScalarImp = double>
class EigenDenseMatrix : public MatrixInterface<EigenDenseMatrixTraits<ScalarImp>>,
                         public EigenMatrixInterfaceDynamic,
                         public ProvidesBackend<EigenDenseMatrixTraits<ScalarImp>>,
                         public ProvidesDataAccess<EigenDenseMatrixTraits<ScalarImp>>
{
  typedef EigenDenseMatrix<ScalarImp> ThisType;
  typedef MatrixInterface<EigenDenseMatrixTraits<ScalarImp>> MatrixInterfaceType;
  static_assert(!std::is_same<DUNE_STUFF_SSIZE_T, int>::value,
                "You have to manually disable the constructor below which uses DUNE_STUFF_SSIZE_T!");

public:
  typedef EigenDenseMatrixTraits<ScalarImp> Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::ScalarType ScalarType;

  EigenDenseMatrix(const size_t rr = 0, const size_t cc = 0, const ScalarType value = ScalarType(0))
    : backend_(new BackendType(rr, cc))
  {
    if (FloatCmp::eq(value, ScalarType(0)))
      backend_->setZero();
    else {
      backend_->setOnes();
      backend_->operator*=(value);
    }
  }

  /// This constructor is needed for the python bindings.
  EigenDenseMatrix(const DUNE_STUFF_SSIZE_T rr, const DUNE_STUFF_SSIZE_T cc = 0, const ScalarType value = ScalarType(0))
    : backend_(new BackendType(MatrixInterfaceType::assert_is_size_t_compatible_and_convert(rr),
                               MatrixInterfaceType::assert_is_size_t_compatible_and_convert(cc)))
  {
    if (FloatCmp::eq(value, ScalarType(0)))
      backend_->setZero();
    else {
      backend_->setOnes();
      backend_->operator*=(value);
    }
  }

  EigenDenseMatrix(const int rr, const int cc = 0, const ScalarType value = ScalarType(0))
    : backend_(new BackendType(MatrixInterfaceType::assert_is_size_t_compatible_and_convert(rr),
                               MatrixInterfaceType::assert_is_size_t_compatible_and_convert(cc)))
  {
    if (FloatCmp::eq(value, ScalarType(0)))
      backend_->setZero();
    else {
      backend_->setOnes();
      backend_->operator*=(value);
    }
  }

  /// This constructors ignores the given pattern and initializes the matrix with 0.
  EigenDenseMatrix(const size_t rr, const size_t cc, const SparsityPatternDefault& /*pattern*/)
    : backend_(new BackendType(MatrixInterfaceType::assert_is_size_t_compatible_and_convert(rr),
                               MatrixInterfaceType::assert_is_size_t_compatible_and_convert(cc)))
  {
    backend_->setZero();
  }

  EigenDenseMatrix(const ThisType& other)
    : backend_(other.backend_)
  {
  }

  EigenDenseMatrix(const BackendType& other)
    : backend_(new BackendType(other))
  {
  }

  /**
   *  \note Takes ownership of backend_ptr in the sense that you must not delete it afterwards!
   */
  EigenDenseMatrix(BackendType* backend_ptr)
    : backend_(backend_ptr)
  {
  }

  EigenDenseMatrix(std::shared_ptr<BackendType> backend_ptr)
    : backend_(backend_ptr)
  {
  }

  ThisType& operator=(const ThisType& other)
  {
    backend_ = other.backend_;
    return *this;
  } // ... operator=(...)

  /**
   *  \note Does a deep copy.
   */
  ThisType& operator=(const BackendType& other)
  {
    backend_ = std::make_shared<BackendType>(other);
    return *this;
  } // ... operator=(...)

  /**
   * \defgroup backend ´´These methods are required by the ProvidesBackend interface.``
   * \{
   */

  BackendType& backend()
  {
    ensure_uniqueness();
    return *backend_;
  } // ... backend(...)

  const BackendType& backend() const
  {
    const_cast<ThisType&>(*this).ensure_uniqueness();
    return *backend_;
  } // ... backend(...)
  /**
   * \}
   */

  /**
   * \defgroup data ´´These methods are required by the ProvidesDataAccess interface.``
   * \{
   */

  ScalarType* data()
  {
    return backend_->data();
  }

  /**
   * \}
   */

  /**
   * \defgroup container ´´These methods are required by ContainerInterface.``
   * \{
   */

  ThisType copy() const
  {
    return ThisType(*backend_);
  }

  void scal(const ScalarType& alpha)
  {
    ensure_uniqueness();
    backend_->operator*=(alpha);
  } // ... scal(...)

  void axpy(const ScalarType& alpha, const ThisType& xx)
  {
    if (!has_equal_shape(xx))
      DUNE_THROW_COLORFULLY(
          Exceptions::shapes_do_not_match,
          "The shape of xx (" << xx.rows() << "x" << xx.cols() << ") does not match the shape of this (" << rows()
                              << "x"
                              << cols()
                              << ")!");
    ensure_uniqueness();
    auto& this_ref     = *backend_;
    const auto& xx_ref = *(xx.backend_);
    this_ref += alpha * xx_ref;
  } // ... axpy(...)

  bool has_equal_shape(const ThisType& other) const
  {
    return (rows() == other.rows()) && (cols() == other.cols());
  }
  /**
   * \}
   */

  /**
   * \defgroup matrix_required ´´These methods are required by MatrixInterface.``
   * \{
   */

  inline size_t rows() const
  {
    return backend_->rows();
  }

  inline size_t cols() const
  {
    return backend_->cols();
  }

  template <class SourceType, class RangeType>
  inline void mv(const SourceType& /*xx*/, RangeType& /*yy*/) const
  {
    static_assert(Dune::AlwaysFalse<SourceType>::value, "Not available for this combination of xx and yy!");
  }

  inline void mv(const EigenDenseVector<ScalarType>& xx, EigenDenseVector<ScalarType>& yy) const
  {
    yy.backend_->transpose() = backend_->operator*(*(xx.backend_));
  }

  inline void mv(const EigenDenseVector<ScalarType>& xx, EigenMappedDenseVector<ScalarType>& yy) const
  {
    yy.backend_->transpose() = backend_->operator*(*(xx.backend_));
  }

  inline void mv(const EigenMappedDenseVector<ScalarType>& xx, EigenDenseVector<ScalarType>& yy) const
  {
    yy.backend_->transpose() = backend_->operator*(*(xx.backend_));
  }

  inline void mv(const EigenMappedDenseVector<ScalarType>& xx, EigenMappedDenseVector<ScalarType>& yy) const
  {
    yy.backend_->transpose() = backend_->operator*(*(xx.backend_));
  }

  void add_to_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    assert(ii < rows());
    assert(jj < cols());
    ensure_uniqueness();
    backend_->operator()(ii, jj) += value;
  } // ... add_to_entry(...)

  void set_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    assert(ii < rows());
    assert(jj < cols());
    ensure_uniqueness();
    backend_->operator()(ii, jj) = value;
  } // ... set_entry(...)

  ScalarType get_entry(const size_t ii, const size_t jj) const
  {
    assert(ii < rows());
    assert(jj < cols());
    return backend_->operator()(ii, jj);
  } // ... get_entry(...)

  void clear_row(const size_t ii)
  {
    if (ii >= rows())
      DUNE_THROW_COLORFULLY(Exceptions::index_out_of_range,
                            "Given ii (" << ii << ") is larger than the rows of this (" << rows() << ")!");
    for (size_t jj = 0; jj < cols(); ++jj)
      backend_->operator()(ii, jj) = ScalarType(0);
  } // ... clear_row(...)

  void clear_col(const size_t jj)
  {
    if (jj >= cols())
      DUNE_THROW_COLORFULLY(Exceptions::index_out_of_range,
                            "Given jj (" << jj << ") is larger than the cols of this (" << cols() << ")!");
    for (size_t ii = 0; ii < rows(); ++ii)
      backend_->operator()(ii, jj) = ScalarType(0);
  } // ... clear_col(...)

  void unit_row(const size_t ii)
  {
    if (ii >= rows())
      DUNE_THROW_COLORFULLY(Exceptions::index_out_of_range,
                            "Given ii (" << ii << ") is larger than the rows of this (" << rows() << ")!");
    for (size_t jj = 0; jj < cols(); ++jj)
      backend_->operator()(ii, jj) = ScalarType(0);
    backend_->operator()(ii, ii) = ScalarType(1);
  } // ... unit_row(...)

  void unit_col(const size_t jj)
  {
    if (jj >= cols())
      DUNE_THROW_COLORFULLY(Exceptions::index_out_of_range,
                            "Given jj (" << jj << ") is larger than the cols of this (" << cols() << ")!");
    for (size_t ii = 0; ii < rows(); ++ii)
      backend_->operator()(ii, jj) = ScalarType(0);
    backend_->operator()(jj, jj) = ScalarType(1);
  } // ... unit_col(...)

  /**
   * \}
   */

private:
  inline void ensure_uniqueness()
  {
    if (!backend_.unique())
      backend_ = std::make_shared<BackendType>(*backend_);
  } // ... ensure_uniqueness(...)

  friend class Dune::Pymor::Operators::EigenRowMajorSparseInverse<ScalarType>;
  friend class Dune::Pymor::Operators::EigenRowMajorSparse<ScalarType>;

  std::shared_ptr<BackendType> backend_;
}; // class EigenDenseMatrix


/**
 * \brief Traits for EigenRowMajorSparseMatrix.
 */
template <class ScalarImp = double>
class EigenRowMajorSparseMatrixTraits
{
public:
  typedef ScalarImp ScalarType;
  typedef EigenRowMajorSparseMatrix<ScalarType> derived_type;
  typedef typename ::Eigen::SparseMatrix<ScalarType, ::Eigen::RowMajor> BackendType;
}; // class RowMajorSparseMatrixTraits


/**
 * \brief A sparse matrix implementation of the MatrixInterface with row major memory layout.
 */
template <class ScalarImp = double>
class EigenRowMajorSparseMatrix : public MatrixInterface<EigenRowMajorSparseMatrixTraits<ScalarImp>>,
                                  public EigenMatrixInterfaceDynamic,
                                  public ProvidesBackend<EigenRowMajorSparseMatrixTraits<ScalarImp>>
{
  typedef EigenRowMajorSparseMatrix<ScalarImp> ThisType;
  typedef MatrixInterface<EigenRowMajorSparseMatrixTraits<ScalarImp>> MatrixInterfaceType;
  static_assert(!std::is_same<DUNE_STUFF_SSIZE_T, int>::value,
                "You have to manually disable the constructor below which uses DUNE_STUFF_SSIZE_T!");

public:
  typedef EigenRowMajorSparseMatrixTraits<ScalarImp> Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::ScalarType ScalarType;

  /**
   * \brief This is the constructor of interest which creates a sparse matrix.
   */
  EigenRowMajorSparseMatrix(const size_t rr, const size_t cc, const SparsityPatternDefault& pattern)
    : backend_(new BackendType(rr, cc))
  {
    if (size_t(pattern.size()) != rr)
      DUNE_THROW_COLORFULLY(
          Exceptions::shapes_do_not_match,
          "The size of the pattern (" << pattern.size() << ") does not match the number of rows of this (" << rows()
                                      << ")!");
    for (size_t row = 0; row < size_t(pattern.size()); ++row) {
      backend_->startVec(row);
      const auto& columns = pattern.inner(row);
      for (auto& column : columns)
        backend_->insertBackByOuterInner(row, column);
      // create diagonal entry (insertBackByOuterInner() can not handle empty rows)
      if (columns.size() == 0)
        backend_->insertBackByOuterInner(row, row);
    }
    backend_->finalize();
    backend_->makeCompressed();
  }

  EigenRowMajorSparseMatrix(const size_t rr = 0, const size_t cc = 0)
    : backend_(new BackendType(rr, cc))
  {
  }

  /// This constructor is needed for the python bindings.
  EigenRowMajorSparseMatrix(const DUNE_STUFF_SSIZE_T rr, const DUNE_STUFF_SSIZE_T cc = 0)
    : backend_(new BackendType(MatrixInterfaceType::assert_is_size_t_compatible_and_convert(rr),
                               MatrixInterfaceType::assert_is_size_t_compatible_and_convert(cc)))
  {
  }

  EigenRowMajorSparseMatrix(const int rr, const int cc = 0)
    : backend_(new BackendType(MatrixInterfaceType::assert_is_size_t_compatible_and_convert(rr),
                               MatrixInterfaceType::assert_is_size_t_compatible_and_convert(cc)))
  {
  }

  EigenRowMajorSparseMatrix(const ThisType& other)
    : backend_(other.backend_)
  {
  }

  EigenRowMajorSparseMatrix(const BackendType& other)
    : backend_(new BackendType(other))
  {
  }

  /**
   *  \note Takes ownership of backend_ptr in the sense that you must not delete it afterwards!
   */
  EigenRowMajorSparseMatrix(BackendType* backend_ptr)
    : backend_(backend_ptr)
  {
  }

  EigenRowMajorSparseMatrix(std::shared_ptr<BackendType> backend_ptr)
    : backend_(backend_ptr)
  {
  }

  ThisType& operator=(const ThisType& other)
  {
    backend_ = other.backend_;
    return *this;
  } // ... operator=(...)

  /**
   *  \note Does a deep copy.
   */
  ThisType& operator=(const BackendType& other)
  {
    backend_ = std::make_shared<BackendType>(other);
    return *this;
  } // ... operator=(...)

  /**
   * \defgroup backend ´´These methods are required by the ProvidesBackend interface.``
   * \{
   */

  BackendType& backend()
  {
    ensure_uniqueness();
    return *backend_;
  } // ... backend(...)

  const BackendType& backend() const
  {
    const_cast<ThisType&>(*this).ensure_uniqueness();
    return *backend_;
  } // ... backend(...)
  /**
   * \}
   */

  /**
   * \defgroup container ´´These methods are required by ContainerInterface.``
   * \{
   */

  ThisType copy() const
  {
    return ThisType(*backend_);
  }

  void scal(const ScalarType& alpha)
  {
    ensure_uniqueness();
    backend_->operator*=(alpha);
  } // ... scal(...)

  void axpy(const ScalarType& alpha, const ThisType& xx)
  {
    if (!has_equal_shape(xx))
      DUNE_THROW_COLORFULLY(
          Exceptions::shapes_do_not_match,
          "The shape of xx (" << xx.rows() << "x" << xx.cols() << ") does not match the shape of this (" << rows()
                              << "x"
                              << cols()
                              << ")!");
    ensure_uniqueness();
    auto& this_ref     = *backend_;
    const auto& xx_ref = *(xx.backend_);
    this_ref += alpha * xx_ref;
  } // ... axpy(...)

  bool has_equal_shape(const ThisType& other) const
  {
    return (rows() == other.rows()) && (cols() == other.cols());
  }
  /**
   * \}
   */

  /**
   * \defgroup matrix_required ´´These methods are required by MatrixInterface.``
   * \{
   */

  inline size_t rows() const
  {
    return backend_->rows();
  }

  inline size_t cols() const
  {
    return backend_->cols();
  }

  //  template< class SourceType, class RangeType >
  //  inline void mv(const SourceType& /*xx*/, RangeType& /*yy*/) const
  //  {
  //    static_assert(Dune::AlwaysFalse< SourceType >::value, "Not available for this combination of xx and yy!");
  //  }

  template <class T1, class T2>
  inline void mv(const EigenBaseVector<T1, ScalarType>& xx, EigenBaseVector<T2, ScalarType>& yy) const
  {
    yy.backend().transpose() = backend_->operator*(xx.backend());
  }

  //  inline void mv(const EigenDenseVector< ScalarType >& xx, EigenMappedDenseVector< ScalarType >& yy) const
  //  {
  //    yy.backend_->transpose() = backend_->operator*(*(xx.backend_));
  //  }

  //  inline void mv(const EigenMappedDenseVector< ScalarType >& xx, EigenDenseVector< ScalarType >& yy) const
  //  {
  //    yy.backend_->transpose() = backend_->operator*(*(xx.backend_));
  //  }

  //  inline void mv(const EigenMappedDenseVector< ScalarType >& xx, EigenMappedDenseVector< ScalarType >& yy) const
  //  {
  //    yy.backend_->transpose() = backend_->operator*(*(xx.backend_));
  //  }

  void add_to_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    assert(these_are_valid_indices(ii, jj));
    backend_->coeffRef(ii, jj) += value;
  } // ... add_to_entry(...)

  void set_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    assert(these_are_valid_indices(ii, jj));
    backend_->coeffRef(ii, jj) = value;
  } // ... set_entry(...)

  ScalarType get_entry(const size_t ii, const size_t jj) const
  {
    assert(ii < rows());
    assert(jj < cols());
    return backend_->coeff(ii, jj);
  } // ... get_entry(...)

  void clear_row(const size_t ii)
  {
    if (ii >= rows())
      DUNE_THROW_COLORFULLY(Exceptions::index_out_of_range,
                            "Given ii (" << ii << ") is larger than the rows of this (" << rows() << ")!");
    backend_->row(ii) *= ScalarType(0);
  } // ... clear_row(...)

  void clear_col(const size_t jj)
  {
    if (jj >= cols())
      DUNE_THROW_COLORFULLY(Exceptions::index_out_of_range,
                            "Given jj (" << jj << ") is larger than the cols of this (" << cols() << ")!");
    for (size_t row = 0; row < backend_->outerSize(); ++row) {
      for (typename BackendType::InnerIterator row_it(*backend_, row); row_it; ++row_it) {
        const size_t col = row_it.col();
        if (col == jj) {
          backend_->coeffRef(row, jj) = ScalarType(0);
          break;
        } else if (col > jj)
          break;
      }
    }
  } // ... clear_col(...)

  void unit_row(const size_t ii)
  {
    if (ii >= rows())
      DUNE_THROW_COLORFULLY(Exceptions::index_out_of_range,
                            "Given ii (" << ii << ") is larger than the rows of this (" << rows() << ")!");
    if (!these_are_valid_indices(ii, ii))
      DUNE_THROW_COLORFULLY(Exceptions::index_out_of_range,
                            "Diagonal entry (" << ii << ", " << ii << ") is not contained in the sparsity pattern!");
    backend_->row(ii) *= ScalarType(0);
    set_entry(ii, ii, ScalarType(1));
  } // ... unit_row(...)

  void unit_col(const size_t jj)
  {
    if (jj >= cols())
      DUNE_THROW_COLORFULLY(Exceptions::index_out_of_range,
                            "Given jj (" << jj << ") is larger than the cols of this (" << cols() << ")!");
    for (size_t row = 0; row < backend_->outerSize(); ++row) {
      for (typename BackendType::InnerIterator row_it(*backend_, row); row_it; ++row_it) {
        const size_t col = row_it.col();
        if (col == jj) {
          if (col == row)
            backend_->coeffRef(row, col) = ScalarType(1);
          else
            backend_->coeffRef(row, jj) = ScalarType(0);
          break;
        } else if (col > jj)
          break;
      }
    }
  } // ... unit_col(...)
  /**
   * \}
   */

private:
  bool these_are_valid_indices(const size_t ii, const size_t jj)
  {
    if (ii >= rows())
      return false;
    if (jj >= cols())
      return false;
    for (size_t row = ii; row < backend_->outerSize(); ++row) {
      for (typename BackendType::InnerIterator row_it(*backend_, row); row_it; ++row_it) {
        const size_t col = row_it.col();
        if ((ii == row) && (jj == col))
          return true;
        else if ((row > ii) && (col > jj))
          return false;
      }
    }
    return false;
  } // ... these_are_valid_indices(...)

  inline void ensure_uniqueness()
  {
    if (!backend_.unique())
      backend_ = std::make_shared<BackendType>(*backend_);
  } // ... ensure_uniqueness(...)

  friend class Dune::Pymor::Operators::EigenRowMajorSparseInverse<ScalarType>;
  friend class Dune::Pymor::Operators::EigenRowMajorSparse<ScalarType>;

  std::shared_ptr<BackendType> backend_;
}; // class EigenRowMajorSparseMatrix


#else // HAVE_EIGEN


template <class VectorImpTraits, class ScalarImp>
class EigenBaseVector
{
  static_assert(Dune::AlwaysFalse<ScalarImp>::value, "You are missing Eigen!");
};

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

template <class ScalarImp>
class EigenRowMajorSparseMatrix
{
  static_assert(Dune::AlwaysFalse<ScalarImp>::value, "You are missing Eigen!");
};


#endif // HAVE_EIGEN

} // namespace LA
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_LA_CONTAINER_EIGEN_HH
