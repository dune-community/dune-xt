// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2014 - 2018)
//   Tobias Leibner  (2014 - 2015, 2017 - 2018)

#ifndef DUNE_XT_LA_CONTAINER_CONTAINER_INTERFACE_HH
#define DUNE_XT_LA_CONTAINER_CONTAINER_INTERFACE_HH

#include <cmath>
#include <limits>
#include <type_traits>

#include <boost/numeric/conversion/cast.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/lock_algorithms.hpp>

#include <dune/xt/common/crtp.hh>
#include <dune/xt/common/math.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/type_traits.hh>

namespace Dune {
namespace XT {
namespace LA {


enum class Backends
{
  common_dense,
  common_sparse,
  istl_dense,
  istl_sparse,
  eigen_dense,
  eigen_sparse,
  none
}; // enum class Backends

static constexpr Backends default_backend =
#if HAVE_EIGEN
    Backends::eigen_sparse;
#else
    Backends::istl_sparse;
#endif

static constexpr Backends default_sparse_backend =
#if HAVE_EIGEN
    Backends::eigen_sparse;
#else
    Backends::istl_sparse;
#endif

static constexpr Backends default_dense_backend =
#if HAVE_EIGEN
    Backends::eigen_dense;
#else
    Backends::common_dense;
#endif


namespace internal {


/**
 * \brief Tries a boost::numeric_cast and throws an Exceptions::wrong_input_given on failure.
 *
 *        This can be used in the ctor initializer list.
 */
template <class Out, class In>
static Out boost_numeric_cast(const In& in)
{
  try {
    return boost::numeric_cast<Out>(in);
  } catch (boost::bad_numeric_cast& ee) {
    DUNE_THROW(Common::Exceptions::wrong_input_given,
               "There was an error in boost converting '" << in << "' to '" << Common::Typename<Out>::value() << "': "
                                                          << ee.what());
  }
} // ... boost_numeric_cast(...)

struct VectorLockGuard
{
  VectorLockGuard(std::shared_ptr<std::vector<std::mutex>>& mutexes)
    : mutexes_(mutexes ? mutexes.get() : nullptr)
  {
    if (mutexes_)
      boost::lock(mutexes_->begin(), mutexes_->end());
  }

  ~VectorLockGuard()
  {
    if (mutexes_)
      for (auto& mutex : *mutexes_)
        mutex.unlock();
  }

  std::vector<std::mutex>* mutexes_;
}; // VectorLockGuard

struct LockGuard
{
  LockGuard(std::shared_ptr<std::vector<std::mutex>>& mutexes, const size_t ii)
    : mutexes_(mutexes ? mutexes.get() : nullptr)
    , index_(mutexes_ ? (ii % mutexes_->size()) : 0)
  {
    if (mutexes_)
      mutexes_->operator[](index_).lock();
  }

  ~LockGuard()
  {
    if (mutexes_)
      mutexes_->operator[](index_).unlock();
  }

  std::vector<std::mutex>* mutexes_;
  const size_t index_;
}; // LockGuard


} // namespace internal


template <class Traits>
class ProvidesBackend : public Common::CRTPInterface<ProvidesBackend<Traits>, Traits>
{
public:
  typedef typename Traits::BackendType BackendType;
  static const constexpr Backends backend_type = Traits::backend_type;

  inline BackendType& backend()
  {
    CHECK_CRTP(this->as_imp().backend());
    return this->as_imp().backend();
  }

  inline const BackendType& backend() const
  {
    CHECK_CRTP(this->as_imp().backend());
    return this->as_imp().backend();
  }
}; // class ProvidesBackend

/**
 * \brief Interface for all containers (vectors and matrices).
 *
 * \note  All derived classes are supposed to implement copy-on-write. This can be achieved by internally holding a
 *        shared_prt to the appropriate backend and by passing this shared_prt around on copy, move or assingment. Any
 *        class method that writes to the backend or exposes a reference to the backend is then required to make a deep
 *        copy of the backend, if it is not the sole owner of this resource. This can for instance be achieved by
 *        calling a private method:
\code
  inline void ensure_uniqueness() const
  {
    if (!backend_.unique())
      backend_ = std::make_shared< BackendType >(*backend_);
  }
\endcode
 */
template <class Traits, class ScalarImp = typename Traits::ScalarType>
class ContainerInterface : public Common::CRTPInterface<ContainerInterface<Traits, ScalarImp>, Traits>
{
  typedef Common::CRTPInterface<ContainerInterface<Traits, ScalarImp>, Traits> CRTP;
  static_assert(std::is_same<ScalarImp, typename Traits::ScalarType>::value, "");

public:
  typedef ScalarImp ScalarType;
  typedef typename Traits::RealType RealType;

  using typename CRTP::derived_type;

  virtual ~ContainerInterface()
  {
  }

  /// \name Have to be implemented by a derived class!
  /// \{

  /**
   * \brief   Creates a (deep) copy of the underlying resource
   * \return  A new container
   */
  inline derived_type copy() const
  {
    CHECK_CRTP(this->as_imp().copy());
    return this->as_imp().copy();
  }

  /**
   * \brief BLAS SCAL operation (in-place sclar multiplication).
   * \param alpha The scalar coefficient with which each element of the container is multiplied.
   */
  inline void scal(const ScalarType& alpha)
  {
    CHECK_AND_CALL_CRTP(this->as_imp().scal(alpha));
  }

  /**
   * \brief BLAS AXPY operation.
   * \param alpha The scalar coefficient with which each element of the container is multiplied
   * \param xx    Container that is to be elementwise added.
   */
  inline void axpy(const ScalarType& alpha, const derived_type& xx)
  {
    CHECK_AND_CALL_CRTP(this->as_imp().axpy(alpha, xx));
  }

  /**
   * \brief Test for equal sizes.
   * \param other Container the sizes of which this is to be compared to.
   */
  inline bool has_equal_shape(const derived_type& other) const
  {
    CHECK_CRTP(this->as_imp().has_equal_shape(other));
    return this->as_imp().has_equal_shape(other);
  }

protected:
  inline void ensure_uniqueness()
  {
    CHECK_AND_CALL_CRTP(this->as_imp().ensure_uniqueness());
  }

public:
  /// \}
  /// \name Are provided by the interface for convenience!
  /// \note Those marked as virtual may be implemented more efficiently in a derived class!
  /// \{

  static std::string type_this()
  {
    return Common::Typename<derived_type>::value();
  }

  virtual derived_type& operator*=(const ScalarType& alpha)
  {
    scal(alpha);
    return this->as_imp();
  }

  /**
   *  \brief  Multiplies every component of this by a scalar.
   *  \param  alpha The scalar.
   *  \return The scaled copy of this.
   */
  virtual derived_type operator*(const ScalarType& alpha) const
  {
    auto ret = this->copy();
    ret.scal(alpha);
    return ret;
  }
  /// \}
}; // class ContainerInterface


template <class Traits>
class ProvidesConstContainer : public Common::CRTPInterface<ProvidesConstContainer<Traits>, Traits>
{
public:
  typedef typename Traits::ContainerType ContainerType;

  inline std::shared_ptr<const ContainerType> container() const
  {
    CHECK_CRTP(this->as_imp().container());
    return this->as_imp().container();
  }
}; // class ProvidesConstContainer


template <class Traits>
class ProvidesContainer : public ProvidesConstContainer<Traits>
{
  typedef ProvidesConstContainer<Traits> BaseType;

public:
  typedef typename Traits::ContainerType ContainerType;

  using BaseType::container;

  inline std::shared_ptr<ContainerType> container()
  {
    CHECK_CRTP(this->as_imp().container());
    return this->as_imp().container();
  }
}; // class ProvidesContainer


template <class Traits>
class ProvidesDataAccess : public Common::CRTPInterface<ProvidesDataAccess<Traits>, Traits>
{
public:
  typedef typename Traits::DataType DataType;

  inline DataType* data()
  {
    CHECK_CRTP(this->as_imp().data());
    return this->as_imp().data();
  }

  inline size_t data_size() const
  {
    CHECK_CRTP(this->as_imp().data_size());
    return this->as_imp().data_size();
  }
}; // class ProvidesDataAccess


} // namespace LA
} // namespace XT
} // namespace Dune

#endif // #ifndef DUNE_XT_LA_CONTAINER_CONTAINER_INTERFACE_HH
