// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Tobias Leibner

#ifndef DUNE_STUFF_LA_CONTAINER_CONTAINER_INTERFACE_HH
#define DUNE_STUFF_LA_CONTAINER_CONTAINER_INTERFACE_HH

#ifndef DUNE_STUFF_SSIZE_T
#define DUNE_STUFF_SSIZE_T long int
#endif

#include <cmath>
#include <limits>
#include <type_traits>

#include <dune/stuff/common/crtp.hh>
#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/common/exceptions.hh>

namespace Dune {
namespace Stuff {
namespace LA {


enum class ChooseBackend
{
  common_dense,
  istl_sparse,
  eigen_dense,
  eigen_sparse
}; // enum class ChooseBackend


static constexpr ChooseBackend default_backend =
#if HAVE_EIGEN
    ChooseBackend::eigen_sparse;
#elif HAVE_DUNE_ISTL
    ChooseBackend::istl_sparse;
#else
    ChooseBackend::common_dense;
#endif


static constexpr ChooseBackend default_sparse_backend =
#if HAVE_EIGEN
    ChooseBackend::eigen_sparse;
#elif HAVE_DUNE_ISTL
    ChooseBackend::istl_sparse;
#else
    ChooseBackend::common_dense;
#error "There is no sparse LA backend available!"
#endif


static constexpr ChooseBackend default_dense_backend =
#if HAVE_EIGEN
    ChooseBackend::eigen_dense;
#else
    ChooseBackend::common_dense;
#endif


/**
 *  \brief  Contains tags mostly needed for python bindings.
 */
namespace Tags {


class ContainerInterface
{
};
class ProvidesDataAccess
{
};


} // namespace Tags


template <class Traits>
class ProvidesBackend : public CRTPInterface<ProvidesBackend<Traits>, Traits>
{
public:
  typedef typename Traits::BackendType BackendType;

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
template <class Traits, class ScalarImp>
class ContainerInterface : public Tags::ContainerInterface,
                           public CRTPInterface<ContainerInterface<Traits, ScalarImp>, Traits>
{
  typedef CRTPInterface<ContainerInterface<Traits, ScalarImp>, Traits> CRTP;
  static_assert(std::is_same<ScalarImp, typename Traits::ScalarType>::value, "");

public:
  typedef ScalarImp ScalarType;

  using typename CRTP::derived_type;

  virtual ~ContainerInterface()
  {
  }

  /**
   * \defgroup haveto ´´These methods have to be implemented by a derived class!``
   * \{
   */

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

  /**
   * \}
   */

  /**
   * \defgroup provided ´´These methods are provided by the interface for convenience! Those marked as virtual may be
   * implemented more efficiently in a derived class!``
   * \{
   */

  static std::string type_this()
  {
    return Common::Typename<derived_type>::value();
  }

  virtual derived_type& operator*=(const ScalarType alpha)
  {
    scal(alpha);
    return this->as_imp(*this);
  }

  /**
   * \}
   */

protected:
  template <class SignedSizeType>
  static size_t DUNE_DEPRECATED_MSG("Use boost::numeric_cast instead (24.09.2014)!")
      assert_is_size_t_compatible_and_convert(const SignedSizeType& size)
  {
    if (size < 0)
      DUNE_THROW(Exceptions::index_out_of_range, "Given size (" << size << ") has to be non-negative!");
    typedef typename std::make_unsigned<SignedSizeType>::type UnsignedSizeType;
    if (UnsignedSizeType(size) > std::numeric_limits<size_t>::max())
      DUNE_THROW(Exceptions::index_out_of_range,
                 "Given size (" << size << ") is too large for size_t (max " << std::numeric_limits<size_t>::max()
                                << ")!");
    return size_t(size);
  } // ... ssize_t_is_valid(...)
}; // class ContainerInterface


template <class Traits>
class ProvidesConstContainer : public CRTPInterface<ProvidesConstContainer<Traits>, Traits>
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
class ProvidesDataAccess : public CRTPInterface<ProvidesDataAccess<Traits>, Traits>, public Tags::ProvidesDataAccess
{
public:
  typedef typename Traits::ScalarType ScalarType;

  inline ScalarType* data()
  {
    CHECK_CRTP(this->as_imp().data());
    return this->as_imp().data();
  }
}; // class ProvidesDataAccess


} // namespace LA
} // namespace Stuff
} // namespace Dune

#endif // #ifndef DUNE_STUFF_LA_CONTAINER_CONTAINER_INTERFACE_HH
