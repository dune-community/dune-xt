// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Tobias Leibner

#ifndef DUNE_STUFF_LA_CONTAINER_VECTOR_INTERFACE_HH
#define DUNE_STUFF_LA_CONTAINER_VECTOR_INTERFACE_HH

#ifndef DUNE_STUFF_SSIZE_T
#define DUNE_STUFF_SSIZE_T long int
#endif

#include <cmath>
#include <limits>
#include <iostream>
#include <type_traits>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/stuff/common/disable_warnings.hh>
#include <dune/common/float_cmp.hh>
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/common/crtp.hh>
#include <dune/stuff/common/exceptions.hh>

#include "container-interface.hh"
#include "vector-interface-internal.hh"


namespace Dune {
namespace Stuff {
namespace LA {


/**
 *  \brief  Contains tags mostly needed for python bindings.
 */
namespace Tags {


class VectorInterface
{
};


} // namespace Tags


template <class Traits, class ScalarImp>
class VectorInterface : public ContainerInterface<Traits, ScalarImp>, public Tags::VectorInterface
{
public:
  typedef typename Traits::derived_type derived_type;
  typedef typename Traits::ScalarType TraitsScalarType;
  typedef ScalarImp ScalarType;
  static_assert(std::is_same<ScalarType, TraitsScalarType>::value, "");

  typedef internal::VectorInputIterator<Traits, ScalarType> const_iterator;
  typedef internal::VectorOutputIterator<Traits, ScalarType> iterator;

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
    CHECK_CRTP(this->as_imp().size());
    return this->as_imp().size();
  } // ... size(...)

  /**
   * \brief Add a scalar to the iith entry.
   */
  inline void add_to_entry(const size_t ii, const ScalarType& value)
  {
    CHECK_AND_CALL_CRTP(this->as_imp().add_to_entry(ii, value));
  }

  /**
   * \brief Set the iith entry to given scalar.
   */
  inline void set_entry(const size_t ii, const ScalarType& value)
  {
    CHECK_AND_CALL_CRTP(this->as_imp().set_entry(ii, value));
  }

  /**
   * \brief Get the iith entry.
   */
  inline ScalarType get_entry(const size_t ii) const
  {
    CHECK_CRTP(this->as_imp().get_entry(ii));
    return this->as_imp().get_entry(ii);
  } // ... get_entry(...)

  inline ScalarType& get_entry_ref(const size_t ii)
  {
    CHECK_CRTP(this->as_imp().get_entry_ref(ii));
    return this->as_imp().get_entry_ref(ii);
  } // ... get_entry_ref(...)

  inline const ScalarType& get_entry_ref(const size_t ii) const
  {
    CHECK_CRTP(this->as_imp().get_entry_ref(ii));
    return this->as_imp().get_entry_ref(ii);
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

  virtual ScalarType mean() const
  {
    ScalarType ret = 0.0;
    for (const auto& element : *this)
      ret += element;
    ret /= size();
    return ret;
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
      DUNE_THROW(Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
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
      DUNE_THROW(Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
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

  virtual ScalarType standard_deviation() const
  {
    const ScalarType mu = mean();
    ScalarType sigma = 0.0;
    for (const auto& x_i : *this)
      sigma += std::pow(x_i - mu, 2);
    sigma /= size();
    return std::sqrt(sigma);
  } // ... standard_deviation(...)

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
      DUNE_THROW(Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    if (result.size() != size())
      DUNE_THROW(Exceptions::shapes_do_not_match,
                 "The size of result (" << result.size() << ") does not match the size of this (" << size() << ")!");
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
      DUNE_THROW(Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
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
      DUNE_THROW(Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    if (result.size() != size())
      DUNE_THROW(Exceptions::shapes_do_not_match,
                 "The size of result (" << result.size() << ") does not match the size of this (" << size() << ")!");
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
      DUNE_THROW(Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    for (size_t ii = 0; ii < size(); ++ii)
      set_entry(ii, get_entry(ii) - other.get_entry(ii));
  } // ... isub(...)

  virtual derived_type& operator=(const ScalarType& value)
  {
    for (auto& element : *this)
      element = value;
    return this->as_imp(*this);
  }

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

  virtual derived_type& operator+=(const ScalarType& scalar)
  {
    for (auto& element : *this)
      element += scalar;
    return this->as_imp(*this);
  }

  virtual derived_type& operator-=(const ScalarType& scalar)
  {
    for (auto& element : *this)
      element -= scalar;
    return this->as_imp(*this);
  }

  virtual derived_type& operator/=(const ScalarType& scalar)
  {
    for (auto& element : *this)
      element /= scalar;
    return this->as_imp(*this);
  }

  virtual derived_type& operator*=(const ScalarType& scalar)
  {
    for (auto& element : *this)
      element *= scalar;
    return this->as_imp(*this);
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
    try {
      return boost::numeric_cast<DUNE_STUFF_SSIZE_T>(dim());
    } catch (boost::bad_numeric_cast& ee) {
      DUNE_THROW(Exceptions::external_error,
                 "There was an error in boost converting '" << dim() << "' to '"
                                                            << Common::Typename<ScalarType>::value()
                                                            << "': "
                                                            << ee.what());
    }
  } // ... pb_dim(...)

  /**
   * \brief Variant of add_to_entry() needed for the python bindings.
   * \see   add_to_entry()
   */
  inline void pb_add_to_entry(const DUNE_STUFF_SSIZE_T ii, const ScalarType& value)
  {
    try {
      add_to_entry(boost::numeric_cast<size_t>(ii), value);
    } catch (boost::bad_numeric_cast& ee) {
      DUNE_THROW(Exceptions::external_error,
                 "There was an error in boost converting '" << ii << "' to '" << Common::Typename<size_t>::value()
                                                            << "': "
                                                            << ee.what());
    }
  } // ... pb_add_to_entry(...)

  /**
   * \brief Variant of set_entry() needed for the python bindings.
   * \see   set_entry()
   */
  inline void pb_set_entry(const DUNE_STUFF_SSIZE_T ii, const ScalarType& value)
  {
    try {
      set_entry(boost::numeric_cast<size_t>(ii), value);
    } catch (boost::bad_numeric_cast& ee) {
      DUNE_THROW(Exceptions::external_error,
                 "There was an error in boost converting '" << ii << "' to '" << Common::Typename<size_t>::value()
                                                            << "': "
                                                            << ee.what());
    }
  } // ... pb_set_entry(...)

  /**
   * \brief Variant of get_entry() needed for the python bindings.
   * \see   get_entry()
   */
  inline ScalarType pb_get_entry(const DUNE_STUFF_SSIZE_T ii)
  {
    try {
      return get_entry(boost::numeric_cast<size_t>(ii));
    } catch (boost::bad_numeric_cast& ee) {
      DUNE_THROW(Exceptions::external_error,
                 "There was an error in boost converting '" << ii << "' to '" << Common::Typename<size_t>::value()
                                                            << "': "
                                                            << ee.what());
    }
  } // ... pb_get_entry(...)

  /**
   * \brief Variant of amax() needed for the python bindings.
   * \see   amax()
   */
  std::vector<ScalarType> pb_amax() const
  {
    const auto max = amax();
    try {
      return {boost::numeric_cast<ScalarType>(max.first), max.second};
    } catch (boost::bad_numeric_cast& ee) {
      DUNE_THROW(Exceptions::external_error,
                 "There was an error in boost converting '" << max.first << "' to '"
                                                            << Common::Typename<ScalarType>::value()
                                                            << "': "
                                                            << ee.what());
    }
  } // ... pb_amax(...)

  std::vector<ScalarType> components(const std::vector<DUNE_STUFF_SSIZE_T>& component_indices) const
  {
    if (component_indices.size() > dim())
      DUNE_THROW(Exceptions::index_out_of_range,
                 "size of component_indices (" << component_indices.size() << ") is larger than the dim of this ("
                                               << dim()
                                               << ")!");
    std::vector<ScalarType> values(component_indices.size(), ScalarType(0));
    try {
      for (size_t ii = 0; ii < component_indices.size(); ++ii) {
        const size_t component = boost::numeric_cast<size_t>(component_indices[ii]);
        if (component >= dim())
          DUNE_THROW(Exceptions::index_out_of_range,
                     "component_indices[" << ii << "] is too large for this (" << dim() << ")!");
        values[ii] = get_entry(component);
      }
    } catch (boost::bad_numeric_cast& ee) {
      DUNE_THROW(Exceptions::external_error, "There was an error in boost during a numeric_cast: " << ee.what());
    }
    return values;
  } // components(...)
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
  template <class T, class S>
  friend std::ostream& operator<<(std::ostream& /*out*/, const VectorInterface<T, S>& /*vector*/);
}; // class VectorInterface


template <class T, class S>
std::ostream& operator<<(std::ostream& out, const VectorInterface<T, S>& vector)
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


} // namespace LA
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_LA_CONTAINER_VECTOR_INTERFACE_HH
