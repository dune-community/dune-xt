// This file is part of the dune-stuff project:
//   https://users.dune-project.org/projects/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Tobias Leibner

#ifndef DUNE_STUFF_LA_CONTAINER_INTERFACE_HH
#define DUNE_STUFF_LA_CONTAINER_INTERFACE_HH

#ifndef DUNE_STUFF_SSIZE_T
#define DUNE_STUFF_SSIZE_T long int
#endif

#include <limits>
#include <iostream>
#include <iterator>

#include <dune/common/float_cmp.hh>

#include <dune/stuff/common/crtp.hh>
#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/common/exceptions.hh>

namespace Dune {
namespace Stuff {
namespace LA {


// forward
template <class Traits>
class VectorInterface;


enum class ChooseBackend
{
  common_dense,
  istl_sparse,
  eigen_dense,
  eigen_sparse
}; // enum class ChooseBackend


static constexpr ChooseBackend default_backend =
#if HAVE_DUNE_ISTL
    ChooseBackend::istl_sparse;
#elif HAVE_EIGEN
    ChooseBackend::eigen_sparse;
#else
    ChooseBackend::common_dense;
#endif


static constexpr ChooseBackend default_sparse_backend =
#if HAVE_DUNE_ISTL
    ChooseBackend::istl_sparse;
#elif HAVE_EIGEN
    ChooseBackend::eigen_sparse;
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
class VectorInterface
{
};
class MatrixInterface
{
};


} // namespace Tagss
namespace internal {


template <class Traits>
class VectorInputIterator : public std::iterator<std::input_iterator_tag, typename Traits::ScalarType>
{
  typedef VectorInputIterator<Traits> ThisType;

public:
  typedef VectorInterface<Traits> VectorType;
  typedef typename VectorType::ScalarType ScalarType;

private:
  struct ConstHolder
  {
    ConstHolder(const VectorType& vec)
      : element(vec)
    {
    }

    const VectorType& element;
  };

public:
  VectorInputIterator(const VectorType& vec, const bool end = false)
    : const_holder_(std::make_shared<ConstHolder>(vec))
    , position_(0)
    , end_(end)
  {
  }

  ThisType& operator++()
  {
    if (!end_ && position_ < (const_holder_->element.size() - 1))
      ++position_;
    else
      end_ = true;
    return *this;
  }

  bool operator==(const ThisType& other)
  {
    return (end_ && other.end_) || ((!end_ && !other.end_) && (position_ == other.position_));
  }

  bool operator!=(const ThisType& other)
  {
    return !operator==(other);
  }

  const ScalarType& operator*() const
  {
    if (end_)
      DUNE_THROW_COLORFULLY(Exceptions::you_are_using_this_wrong, "This is the end!");
    return const_holder_->element[position_];
  }

private:
  std::shared_ptr<ConstHolder> const_holder_;

protected:
  size_t position_;
  bool end_;
}; // class VectorInputIterator


template <class Traits>
class VectorOutputIterator : public VectorInputIterator<Traits>,
                             public std::iterator<std::output_iterator_tag, typename Traits::ScalarType>
{
  typedef VectorInputIterator<Traits> BaseType;
  typedef VectorOutputIterator<Traits> ThisType;

public:
  typedef VectorInterface<Traits> VectorType;
  typedef typename VectorType::ScalarType ScalarType;

private:
  struct Holder
  {
    Holder(VectorType& vec)
      : element(vec)
    {
    }

    VectorType& element;
  };

public:
  VectorOutputIterator(VectorType& vec, const bool end = false)
    : BaseType(vec, end)
    , holder_(std::make_shared<Holder>(vec))
  {
  }

  ScalarType& operator*()
  {
    if (this->end_)
      DUNE_THROW_COLORFULLY(Exceptions::you_are_using_this_wrong, "This is the end!");
    return holder_->element[this->position_];
  }

private:
  std::shared_ptr<Holder> holder_;
}; // class VectorOutputIterator


} // namespace internal


template <class Traits>
class ProvidesBackend : public CRTPInterface<ProvidesBackend<Traits>, Traits>
{
public:
  typedef typename Traits::BackendType BackendType;

  inline BackendType& backend()
  {
    CHECK_CRTP(this->as_imp(*this).backend());
    return this->as_imp(*this).backend();
  }

  inline const BackendType& backend() const
  {
    CHECK_CRTP(this->as_imp(*this).backend());
    return this->as_imp(*this).backend();
  }
}; // class ProvidesBackend


template <class Traits>
class ContainerInterface : public Tags::ContainerInterface, public CRTPInterface<ContainerInterface<Traits>, Traits>
{
  typedef CRTPInterface<ContainerInterface<Traits>, Traits> CRTP;

public:
  typedef typename Traits::ScalarType ScalarType;

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
    CHECK_CRTP(this->as_imp(*this).copy());
    return this->as_imp(*this).copy();
  }

  /**
   * \brief BLAS SCAL operation (in-place sclar multiplication).
   * \param alpha The scalar coefficient with which each element of the container is multiplied.
   */
  inline void scal(const ScalarType& alpha)
  {
    CHECK_AND_CALL_CRTP(this->as_imp(*this).scal(alpha));
  }

  /**
   * \brief BLAS AXPY operation.
   * \param alpha The scalar coefficient with which each element of the container is multiplied
   * \param xx    Container that is to be elementwise added.
   */
  inline void axpy(const ScalarType& alpha, const derived_type& xx)
  {
    CHECK_AND_CALL_CRTP(this->as_imp(*this).axpy(alpha, xx));
  }

  /**
   * \brief Test for equal sizes.
   * \param other Container the sizes of which this is to be compared to.
   */
  inline bool has_equal_shape(const derived_type& other) const
  {
    CHECK_CRTP(this->as_imp(*this).has_equal_shape(other));
    return this->as_imp(*this).has_equal_shape(other);
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
    return Stuff::Common::Typename<derived_type>::value();
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
  static size_t assert_is_size_t_compatible_and_convert(const SignedSizeType& size)
  {
    if (size < 0)
      DUNE_THROW_COLORFULLY(Exceptions::index_out_of_range, "Given size (" << size << ") has to be non-negative!");
    if (size > std::numeric_limits<size_t>::max())
      DUNE_THROW_COLORFULLY(
          Exceptions::index_out_of_range,
          "Given size (" << size << ") is to large for size_t (max " << std::numeric_limits<size_t>::max() << ")!");
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
    CHECK_CRTP(this->as_imp(*this).container());
    return this->as_imp(*this).container();
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
    CHECK_CRTP(this->as_imp(*this).container());
    return this->as_imp(*this).container();
  }
}; // class ProvidesContainer


template <class Traits>
class ProvidesDataAccess : public CRTPInterface<ProvidesDataAccess<Traits>, Traits>, public Tags::ProvidesDataAccess
{
public:
  typedef typename Traits::ScalarType ScalarType;

  inline ScalarType* data()
  {
    CHECK_CRTP(this->as_imp(*this).data());
    return this->as_imp(*this).data();
  }
}; // class ProvidesDataAccess


template <class Traits>
class VectorInterface : public ContainerInterface<Traits>, public Tags::VectorInterface
{
public:
  typedef typename Traits::derived_type derived_type;
  typedef typename Traits::ScalarType ScalarType;

  typedef internal::VectorInputIterator<Traits> const_iterator;
  typedef internal::VectorOutputIterator<Traits> iterator;

  virtual ~VectorInterface()
  {
  }

  /**
   * \defgroup haveto ´´These methods have to be implemented by a derived class in addition to the ones required by
   * ContainerInterface!``
   * \{
   */

  /**
   * \brief   The size of the vector.
   * \return  The size of the vector.
   */
  inline size_t size() const
  {
    CHECK_CRTP(this->as_imp(*this).size());
    return this->as_imp(*this).size();
  } // ... size(...)

  /**
   * \brief Add a scalar to the iith entry.
   */
  inline void add_to_entry(const size_t ii, const ScalarType& value)
  {
    CHECK_AND_CALL_CRTP(this->as_imp(*this).add_to_entry(ii, value));
  }

  /**
   * \brief Set the iith entry to given scalar.
   */
  inline void set_entry(const size_t ii, const ScalarType& value)
  {
    CHECK_AND_CALL_CRTP(this->as_imp(*this).set_entry(ii, value));
  }

  /**
   * \brief Get the iith entry.
   */
  inline ScalarType get_entry(const size_t ii) const
  {
    CHECK_CRTP(this->as_imp(*this).get_entry(ii));
    return this->as_imp(*this).get_entry(ii);
  } // ... get_entry(...)

  inline ScalarType& get_entry_ref(const size_t ii)
  {
    CHECK_CRTP(this->as_imp(*this).get_entry_ref(ii));
    return this->as_imp(*this).get_entry_ref(ii);
  } // ... get_entry_ref(...)

  inline const ScalarType& get_entry_ref(const size_t ii) const
  {
    CHECK_CRTP(this->as_imp(*this).get_entry_ref(ii));
    return this->as_imp(*this).get_entry_ref(ii);
  } // ... get_entry_ref(...)

  /**
   * \}
   */

  /**
   * \defgroup provided ´´These methods are provided by the interface for convenience! Those marked as virtual may be
   * implemented more efficiently in a derived class!``
   * \{
   */

  virtual bool valid() const
  {
    for (const auto& val : *this) {
      if (std::isnan(val) || std::isinf(val))
        return false;
    }
    return true;
  }

  /**
   * \brief Get reference to the iith entry.
   */
  inline ScalarType& operator[](const size_t ii)
  {
    return get_entry_ref(ii);
  }

  inline const ScalarType& operator[](const size_t ii) const
  {
    return get_entry_ref(ii);
  }

  /**
   * \brief   The dimension of the vector.
   * \return  The dimension of the vector.
   * \see     size()
   */
  inline size_t dim() const
  {
    return size();
  }

  /**
   *  \brief  The maximum absolute value of the vector.
   *  \return A pair of the lowest index at which the maximum is attained and the absolute maximum value.
   *  \note   If you override this method please use exceptions instead of assertions (for the python bindings).
   */
  virtual std::pair<size_t, ScalarType> amax() const
  {
    auto result = std::make_pair(size_t(0), ScalarType(0));
    for (size_t ii = 0; ii < size(); ++ii) {
      const auto value = std::abs(get_entry(ii));
      if (value > result.second) {
        result.first  = ii;
        result.second = value;
      }
    }
    return result;
  } // ... amax(...)

  /**
   *  \brief  Check vectors for equality.
   *          Equality of two vectors is defined as in Dune::FloatCmp componentwise.
   *  \param  other   A vector of same dimension to compare with.
   *  \param  epsilon See Dune::FloatCmp.
   *  \return Truth value of the comparison.
   *  \see    Dune::FloatCmp
   *  \note   If you override this method please use exceptions instead of assertions (for the python bindings).
   */
  virtual bool almost_equal(const derived_type& other,
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

  /**
   *  \brief  Computes the scalar products between two vectors.
   *  \param  other The second factor.
   *  \return The scalar product.
   *  \note   If you override this method please use exceptions instead of assertions (for the python bindings).
   */
  virtual ScalarType dot(const derived_type& other) const
  {
    if (other.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of other (" << other.size() << ") does not match the size of this (" << size()
                                                  << ")!");
    ScalarType result = 0;
    for (size_t ii = 0; ii < size(); ++ii)
      result += get_entry(ii) * other.get_entry(ii);
    return result;
  } // ... dot(...)

  /**
   *  \brief  The l1-norm of the vector.
   *  \return The l1-norm of the vector.
   *  \note   If you override this method please use exceptions instead of assertions (for the python bindings).
   */
  virtual ScalarType l1_norm() const
  {
    ScalarType result = 0;
    for (size_t ii = 0; ii < size(); ++ii)
      result += std::abs(get_entry(ii));
    return result;
  } // ... l1_norm(...)

  /**
   *  \brief  The l2-norm of the vector.
   *  \return The l2-norm of the vector.
   *  \note   If you override this method please use exceptions instead of assertions (for the python bindings).
   */
  virtual ScalarType l2_norm() const
  {
    return std::sqrt(dot(this->as_imp(*this)));
  }

  /**
   *  \brief  The l-infintiy-norm of the vector.
   *  \return The l-infintiy-norm of the vector.
   *  \note   If you override this method please use exceptions instead of assertions (for the python bindings).
   */
  virtual ScalarType sup_norm() const
  {
    return amax().second;
  }

  /**
   *  \brief  Adds two vectors.
   *  \param  other   The right summand.
   *  \param  result  Vector to write the result of this + other to
   *  \note   If you override this method please use exceptions instead of assertions (for the python bindings).
   *  \note   If you are looking for the old (now deprecated) add() method, \see add_to_entry().
   */
  virtual void add(const derived_type& other, derived_type& result) const
  {
    if (other.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of other (" << other.size() << ") does not match the size of this (" << size()
                                                  << ")!");
    if (result.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of result (" << result.size() << ") does not match the size of this (" << size()
                                                   << ")!");
    for (size_t ii = 0; ii < size(); ++ii)
      result.set_entry(ii, get_entry(ii) + other.get_entry(ii));
  } // ... add(...)

  /**
   *  \brief  Adds two vectors.
   *  \param  other The right summand.
   *  \return The sum of this and other.
   *  \note   If you override this method please use exceptions instead of assertions (for the python bindings).
   *  \note   If you are looking for the old (now deprecated) add() method, \see add_to_entry().
   */
  virtual derived_type add(const derived_type& other) const
  {
    derived_type result = this->copy();
    result.iadd(other);
    return result;
  } // ... add(...)

  /**
   *  \brief  Inplace variant of add().
   *  \param  other The right summand.
   *  \note   If you override this method please use exceptions instead of assertions (for the python bindings).
   */
  virtual void iadd(const derived_type& other)
  {
    if (other.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of other (" << other.size() << ") does not match the size of this (" << size()
                                                  << ")!");
    for (size_t ii = 0; ii < size(); ++ii)
      set_entry(ii, get_entry(ii) + other.get_entry(ii));
  } // ... iadd(...)

  /**
   *  \brief  Subtracts two vectors.
   *  \param  other   The subtrahend.
   *  \param  result  The vectror to write the difference between this and other to.
   *  \note   If you override this method please use exceptions instead of assertions (for the python bindings).
   */
  virtual void sub(const derived_type& other, derived_type& result) const
  {
    if (other.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of other (" << other.size() << ") does not match the size of this (" << size()
                                                  << ")!");
    if (result.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of result (" << result.size() << ") does not match the size of this (" << size()
                                                   << ")!");
    for (size_t ii = 0; ii < size(); ++ii)
      result.set_entry(ii, get_entry(ii) - other.get_entry(ii));
  } // ... sub(...)

  /**
   *  \brief  Subtracts two vectors.
   *  \param  other The subtrahend.
   *  \return The difference between this and other.
   *  \note   If you override this method please use exceptions instead of assertions (for the python bindings).
   */
  virtual derived_type sub(const derived_type& other) const
  {
    derived_type result = this->copy();
    result.isub(other);
    return result;
  } // ... sub(...)

  /**
   *  \brief  Inplace variant of sub().
   *  \param  other The subtrahend.
   *  \note   If you override this method please use exceptions instead of assertions (for the python bindings).
   */
  virtual void isub(const derived_type& other)
  {
    if (other.size() != size())
      DUNE_THROW_COLORFULLY(Exceptions::shapes_do_not_match,
                            "The size of other (" << other.size() << ") does not match the size of this (" << size()
                                                  << ")!");
    for (size_t ii = 0; ii < size(); ++ii)
      set_entry(ii, get_entry(ii) - other.get_entry(ii));
  } // ... isub(...)

  /**
   *  \brief  Multiplies every component of this by a scalar.
   *  \param  alpha The scalar.
   *  \return The scaled copy of this.
   */
  virtual derived_type operator*(const ScalarType& alpha)
  {
    derived_type ret = this->copy();
    ret *= alpha;
    return ret;
  } // ... operator*() ...

  /**
   *  \brief  Computes the scalar products between this and another vector.
   *  \param  other The second factor.
   *  \return The scalar product.
   *  \see dot()
   */
  virtual ScalarType operator*(const derived_type& other)
  {
    return dot(other);
  }

  /**
   *  \brief  Adds another vector to this, in-place variant.
   *  \param  other The second summand.
   *  \return The sum of this and other.
   */
  virtual derived_type& operator+=(const derived_type& other)
  {
    iadd(other);
    return this->as_imp(*this);
  }

  /**
   *  \brief  Subtracts another vector from this, in-place variant.
   *  \param  other The subtrahend.
   *  \return The difference between this and other.
   */
  virtual derived_type& operator-=(const derived_type& other)
  {
    isub(other);
    return this->as_imp(*this);
  }

  /**
   *  \brief  Adds two vectors.
   *  \param  other The second summand.
   *  \return The sum of the two vectors.
   */
  virtual derived_type operator+(const derived_type& other) const
  {
    return add(other);
  }

  /**
   *  \brief  Substracts two vectors.
   *  \param  other The subtrahend
   *  \return The difference.
   */
  virtual derived_type operator-(const derived_type& other) const
  {
    return sub(other);
  }

  /**
   *  \brief  Check vectors for equality (componentwise) using almost_equal()
   *  \param  other   A vector of same dimension to compare with.
   *  \return Truth value of the comparison.
   *  \see    almost_equal()
   */
  virtual bool operator==(const derived_type& other) const
  {
    return almost_equal(other);
  }

  /**
   *  \brief  Check vectors for inequality using !almost_equal()
   *  \param  other   A vector of same dimension to compare with.
   *  \return Truth value of the comparison.
   */
  virtual bool operator!=(const derived_type& other) const
  {
    return !(this->operator==(other));
  }
  /**
   * \}
   */

  /**
   * \defgroup python_bindings ´´These methods are necesarry for the python bindings. They are provided by the
   * interface!``
   * \{
   */

  /**
   * \brief Variant of dim() needed for the python bindings.
   * \see   dim()
   */
  inline DUNE_STUFF_SSIZE_T pb_dim() const
  {
    if (!(dim() <= std::numeric_limits<DUNE_STUFF_SSIZE_T>::max()))
      DUNE_THROW_COLORFULLY(Exceptions::index_out_of_range,
                            "The dimension of this (" << dim() << ") does not fit into DUNE_STUFF_SSIZE_T (max "
                                                      << std::numeric_limits<DUNE_STUFF_SSIZE_T>::max()
                                                      << ")!");
    return (DUNE_STUFF_SSIZE_T)(dim());
  } // ... pb_dim(...)

  /**
   * \brief Variant of add_to_entry() needed for the python bindings.
   * \see   add_to_entry()
   */
  inline void pb_add_to_entry(const DUNE_STUFF_SSIZE_T ii, const ScalarType& value)
  {
    add_to_entry(this->assert_is_size_t_compatible_and_convert(ii), value);
  } // ... pb_add_to_entry(...)

  /**
   * \brief Variant of set_entry() needed for the python bindings.
   * \see   set_entry()
   */
  inline void pb_set_entry(const DUNE_STUFF_SSIZE_T ii, const ScalarType& value)
  {
    set_entry(this->assert_is_size_t_compatible_and_convert(ii), value);
  } // ... pb_set_entry(...)

  /**
   * \brief Variant of get_entry() needed for the python bindings.
   * \see   get_entry()
   */
  inline ScalarType pb_get_entry(const DUNE_STUFF_SSIZE_T ii)
  {
    return get_entry(this->assert_is_size_t_compatible_and_convert(ii));
  } // ... pb_get_entry(...)

  /**
   * \brief Variant of amax() needed for the python bindings.
   * \see   amax()
   */
  std::vector<ScalarType> pb_amax() const
  {
    const auto max = amax();
    if (!(max.first >= std::numeric_limits<ScalarType>::min() && max.first <= std::numeric_limits<ScalarType>::max()))
      DUNE_THROW_COLORFULLY(Exceptions::index_out_of_range,
                            "The index (" << max.first << ") does not fit into ScalarType (min "
                                          << std::numeric_limits<DUNE_STUFF_SSIZE_T>::min()
                                          << ", max "
                                          << std::numeric_limits<ScalarType>::max()
                                          << ")!");
    std::vector<ScalarType> ret(2);
    ret[0] = ScalarType(max.first);
    ret[1] = max.second;
    return ret;
  } // ... pb_amax(...)

  std::vector<ScalarType> components(const std::vector<DUNE_STUFF_SSIZE_T>& component_indices) const
  {
    if (component_indices.size() > dim())
      DUNE_THROW_COLORFULLY(
          Exceptions::index_out_of_range,
          "size of component_indices (" << component_indices.size() << ") is larger than the dim of this (" << dim()
                                        << ")!");
    std::vector<ScalarType> values(component_indices.size(), ScalarType(0));
    for (size_t ii = 0; ii < component_indices.size(); ++ii) {
      const size_t component = this->assert_is_size_t_compatible_and_convert(component_indices[ii]);
      if (component >= dim())
        DUNE_THROW_COLORFULLY(Exceptions::index_out_of_range,
                              "component_indices[" << ii << "] is too large for this (" << dim() << ")!");
      values[ii] = get_entry(component);
    }
    return values;
  } // components(...)
  /**
   * \}
   */

  /**
   * \defgroup deprecated ´´These methods are provided for backwards compatibility only and will be removed soon!``
   * \{
   */
  inline void DUNE_DEPRECATED_MSG("Please use add_to_entry()!") add(const size_t ii, const ScalarType& value)
  {
    add_to_entry(ii, value);
  }

  inline void DUNE_DEPRECATED_MSG("Please use set_entry()!") set(const size_t ii, const ScalarType& value)
  {
    set_entry(ii, value);
  }

  inline ScalarType DUNE_DEPRECATED_MSG("Please use get_entry()!") get(const size_t ii) const
  {
    return get_entry(ii);
  }
  /**
   * \}
   */

  iterator begin()
  {
    return iterator(*this);
  }

  const_iterator begin() const
  {
    return const_iterator(*this);
  }

  iterator end()
  {
    return iterator(*this, true);
  }

  const_iterator end() const
  {
    return const_iterator(*this, true);
  }

private:
  template <class T>
  friend std::ostream& operator<<(std::ostream& /*out*/, const VectorInterface<T>& /*vector*/);
}; // class VectorInterface


template <class T>
std::ostream& operator<<(std::ostream& out, const VectorInterface<T>& vector)
{
  out << "[";
  const size_t sz = vector.size();
  if (sz > 0) {
    out << vector[0];
    for (size_t ii = 1; ii < sz; ++ii)
      out << "\n " << vector[ii];
  } else
    out << " ";
  out << "]";
  return out;
} // ... operator<<(...)


template <class Traits>
class MatrixInterface : public ContainerInterface<Traits>, public Tags::MatrixInterface
{
public:
  typedef typename Traits::derived_type derived_type;
  typedef typename Traits::ScalarType ScalarType;

  virtual ~MatrixInterface()
  {
  }

  /**
   * \defgroup haveto ´´These methods have to be implemented by a derived class in addition to the ones required by
   * ContainerInterface!``
   * \{
   */

  inline size_t rows() const
  {
    CHECK_CRTP(this->as_imp(*this).rows());
    return this->as_imp(*this).rows();
  }

  inline size_t cols() const
  {
    CHECK_CRTP(this->as_imp(*this).cols());
    return this->as_imp(*this).cols();
  }

  template <class XX, class YY>
  inline void mv(const XX& xx, YY& yy) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp(*this).mv(xx, yy));
  }

  inline void add_to_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    CHECK_AND_CALL_CRTP(this->as_imp(*this).add_to_entry(ii, jj, value));
  }

  inline void set_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    CHECK_AND_CALL_CRTP(this->as_imp(*this).set_entry(ii, jj, value));
  }

  inline ScalarType get_entry(const size_t ii, const size_t jj) const
  {
    CHECK_CRTP(this->as_imp(*this).get_entry(ii, jj));
    return this->as_imp(*this).get_entry(ii, jj);
  }

  inline void clear_row(const size_t ii)
  {
    CHECK_AND_CALL_CRTP(this->as_imp(*this).clear_row(ii));
  }

  inline void clear_col(const size_t jj)
  {
    CHECK_AND_CALL_CRTP(this->as_imp(*this).clear_col(jj));
  }

  inline void unit_row(const size_t ii)
  {
    CHECK_AND_CALL_CRTP(this->as_imp(*this).unit_row(ii));
  }

  inline void unit_col(const size_t jj)
  {
    CHECK_AND_CALL_CRTP(this->as_imp(*this).unit_col(jj));
  }

  /**
   * \}
   */

  /**
   * \defgroup python_bindings ´´These methods are necesarry for the python bindings. They are provided by the
   * interface!``
   * \{
   */

  inline DUNE_STUFF_SSIZE_T pb_rows() const
  {
    if (!(rows() <= std::numeric_limits<DUNE_STUFF_SSIZE_T>::max()))
      DUNE_THROW_COLORFULLY(Exceptions::index_out_of_range,
                            "The number of rows of this (" << rows() << ") do not fit into DUNE_STUFF_SSIZE_T (max "
                                                           << std::numeric_limits<DUNE_STUFF_SSIZE_T>::max()
                                                           << ")!");
    return (DUNE_STUFF_SSIZE_T)(rows());
  } // ... pb_rows(...)

  inline DUNE_STUFF_SSIZE_T pb_cols() const
  {
    if (!(cols() <= std::numeric_limits<DUNE_STUFF_SSIZE_T>::max()))
      DUNE_THROW_COLORFULLY(Exceptions::index_out_of_range,
                            "The number of columns of this (" << cols() << ") do not fit into DUNE_STUFF_SSIZE_T (max "
                                                              << std::numeric_limits<DUNE_STUFF_SSIZE_T>::max()
                                                              << ")!");
    return (DUNE_STUFF_SSIZE_T)(cols());
  } // ... pb_cols(...)

  inline void pb_add_to_entry(const DUNE_STUFF_SSIZE_T ii, const DUNE_STUFF_SSIZE_T jj, const ScalarType& value)
  {
    add_to_entry(
        this->assert_is_size_t_compatible_and_convert(ii), this->assert_is_size_t_compatible_and_convert(jj), value);
  } // ... pb_add_to_entry(...)

  inline void pb_set_entry(const DUNE_STUFF_SSIZE_T ii, const DUNE_STUFF_SSIZE_T jj, const ScalarType& value)
  {
    set_entry(
        this->assert_is_size_t_compatible_and_convert(ii), this->assert_is_size_t_compatible_and_convert(jj), value);
  } // ... pb_set_entry(...)

  inline ScalarType pb_get_entry(const DUNE_STUFF_SSIZE_T ii, const DUNE_STUFF_SSIZE_T jj) const
  {
    return get_entry(this->assert_is_size_t_compatible_and_convert(ii),
                     this->assert_is_size_t_compatible_and_convert(jj));
  } // ... pb_get_entry(...)

  inline void pb_clear_row(const size_t ii)
  {
    clear_row(this->assert_is_size_t_compatible_and_convert(ii));
  }

  inline void pb_clear_col(const size_t jj)
  {
    clear_col(this->assert_is_size_t_compatible_and_convert(jj));
  }

  inline void pb_unit_row(const size_t ii)
  {
    unit_row(this->assert_is_size_t_compatible_and_convert(ii));
  }

  inline void pb_unit_col(const size_t jj)
  {
    unit_col(this->assert_is_size_t_compatible_and_convert(jj));
  }

  /**
   * \}
   */

  /**
   * \defgroup deprecated ´´These methods are provided for backwards compatibility only and will be removed soon!``
   * \{
   */
  inline void DUNE_DEPRECATED_MSG("Please use add_to_entry()!")
      add(const size_t ii, const size_t jj, const ScalarType& value)
  {
    add_to_entry(ii, jj, value);
  }

  inline void DUNE_DEPRECATED_MSG("Please use set_entry()!")
      set(const size_t ii, const size_t jj, const ScalarType& value)
  {
    set_entry(ii, jj, value);
  }

  inline ScalarType DUNE_DEPRECATED_MSG("Please use get_entry()!") get(const size_t ii, const size_t jj) const
  {
    return get_entry(ii, jj);
  }

  inline void DUNE_DEPRECATED_MSG("Please use unit_row()!") unitRow(const size_t row)
  {
    unit_row(row);
  }

  inline void DUNE_DEPRECATED_MSG("Please use unit_col()!") unitCol(const size_t col)
  {
    unit_col(col);
  }

  inline void DUNE_DEPRECATED_MSG("Please use unit_col()!") clearRow(const size_t ii)
  {
    clear_row(ii);
  }

  /**
   * \}
   */

private:
  template <class T>
  friend std::ostream& operator<<(std::ostream& /*out*/, const MatrixInterface<T>& /*matrix*/);
}; // class MatrixInterface


template <class T>
std::ostream& operator<<(std::ostream& out, const MatrixInterface<T>& matrix)
{
  out << "[";
  const size_t rows = matrix.rows();
  const size_t cols = matrix.cols();
  if (rows > 0 && cols > 0) {
    for (size_t ii = 0; ii < rows; ++ii) {
      if (ii > 0)
        out << "\n ";
      out << "[" << matrix.get_entry(ii, 0);
      for (size_t jj = 1; jj < cols; ++jj)
        out << " " << matrix.get_entry(ii, jj);
      out << "]";
      if (rows > 1 && ii < (rows - 1))
        out << ",";
    }
    out << "]";
  } else
    out << "[ ]]";
  return out;
} // ... operator<<(...)


} // namespace LA
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_LA_CONTAINER_INTERFACE_HH
