// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Barbara Verfürth (2015)
//   Felix Schindler  (2012 - 2018)
//   Rene Milk        (2013 - 2018)
//   Tobias Leibner   (2014, 2017 - 2018)

#ifndef DUNE_XT_LA_CONTAINER_VECTOR_INTERFACE_HH
#define DUNE_XT_LA_CONTAINER_VECTOR_INTERFACE_HH

#include <cmath>
#include <limits>
#include <iostream>
#include <vector>
#include <complex>

#include <dune/common/bigunsignedint.hh>
#include <dune/common/unused.hh>

#include <dune/xt/common/crtp.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/math.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/common/vector.hh>

#include <dune/xt/la/type_traits.hh>

#include "container-interface.hh"
#include "vector-interface-internal.hh"

namespace Dune {
namespace XT {
namespace LA {


template <class Traits, class ScalarImp = typename Traits::ScalarType>
class VectorInterface : public ContainerInterface<Traits, ScalarImp>
{
  typedef ContainerInterface<Traits, ScalarImp> BaseType;

public:
  using typename BaseType::derived_type;
  typedef typename Dune::FieldTraits<ScalarImp>::field_type ScalarType;
  typedef typename Dune::FieldTraits<ScalarImp>::real_type RealType;
  static const constexpr Backends dense_matrix_type = Traits::dense_matrix_type;
  static const constexpr Backends sparse_matrix_type = Traits::sparse_matrix_type;

  typedef internal::VectorInputIterator<Traits, ScalarType> const_iterator;
  typedef internal::VectorOutputIterator<Traits, ScalarType> iterator;
  friend const_iterator;
  friend iterator;

  static_assert(std::is_same<ScalarType, typename Traits::ScalarType>::value, "");

  virtual ~VectorInterface()
  {
  }

  /// \name Have to be implemented by a derived class in addition to the ones required by ContainerInterface!
  /// \{

  /**
   * \brief   The size of the vector.
   * \return  The size of the vector.
   */
  inline size_t size() const
  {
    CHECK_CRTP(this->as_imp().size());
    return this->as_imp().size();
  }

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
  }

protected:
  /**
   * The purpose of get_unchecked_ref is to allow direct access to the underlying data without any checks (regarding cow
   * or thread safety). This allows default implementations in the interface with locking prior to for-loops.
   */
  inline ScalarType& get_unchecked_ref(const size_t ii)
  {
    CHECK_CRTP(this->as_imp().get_unchecked_ref(ii));
    return this->as_imp().get_unchecked_ref(ii);
  }

  inline const ScalarType& get_unchecked_ref(const size_t ii) const
  {
    CHECK_CRTP(this->as_imp().get_unchecked_ref(ii));
    return this->as_imp().get_unchecked_ref(ii);
  }

public:
  /// \}
  /// \name Provided by the interface for convenience!
  /// \note Those marked as virtual may be implemented more efficiently in a derived class!
  /// \{

  virtual void set_all(const ScalarType& val)
  {
    for (auto& element : *this)
      element = val;
  }

  virtual bool valid() const
  {
    for (const auto& val : *this) {
      if (Common::isnan(val) || Common::isinf(val))
        return false;
    }
    return true;
  } // ... valid()

  /**
   * \brief Get writable reference to the iith entry.
   */
  inline ScalarType& operator[](const size_t ii)
  {
    CHECK_CRTP(this->as_imp()[ii]);
    return this->as_imp()[ii];
  }

  /**
   * \brief Get read-only reference to the iith entry.
   */
  inline const ScalarType& operator[](const size_t ii) const
  {
    CHECK_CRTP(this->as_imp()[ii]);
    return this->as_imp()[ii];
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
  } // ... mean()

  /**
   *  \brief  The maximum absolute value of the vector.
   *  \return A pair of the lowest index at which the maximum is attained and the absolute maximum value.
   *  \note   If you override this method please use exceptions instead of assertions (for the python bindings).
   */
  virtual std::pair<size_t, RealType> amax() const
  {
    using std::abs;
    auto result = std::make_pair(size_t(0), RealType(0));
    for (size_t ii = 0; ii < size(); ++ii) {
      const auto value = abs(get_unchecked_ref(ii));
      if (value > result.second) {
        result.first = ii;
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
   *  \see    Dune::XT::Common::FloatCmp
   *  \note   If you override this method please use exceptions instead of assertions (for the python bindings).
   */
  virtual bool almost_equal(const derived_type& other,
                            const ScalarType epsilon = Common::FloatCmp::DefaultEpsilon<ScalarType>::value()) const
  {
    if (other.size() != size())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    return Common::FloatCmp::eq(this->as_imp(), other, epsilon);
  } // ... almost_equal(...)

  /**
   *  \brief  Check vectors for equality (variant for arbitrary derived combinations).
   *          Equality of two vectors is defined as in Dune::FloatCmp componentwise.
   *  \param  other   A vector of same dimension to compare with.
   *  \param  epsilon See Dune::FloatCmp.
   *  \return Truth value of the comparison.
   *  \see    Dune::XT::Common::FloatCmp
   */
  template <class T>
  bool almost_equal(const VectorInterface<T>& other,
                    const ScalarType epsilon = Common::FloatCmp::DefaultEpsilon<ScalarType>::value()) const
  {
    if (other.size() != size())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    return Common::FloatCmp::eq(this->as_imp(), other.as_imp(), epsilon);
  } // ... almost_equal(...)

private:
  template <bool is_complex = XT::Common::is_complex<ScalarType>::value, bool anything = true>
  struct dot_helper
  {
    static ScalarType compute(const derived_type& ths, const derived_type& other)
    {
      using std::conj;
      if (other.size() != ths.size())
        DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                   "The size of other (" << other.size() << ") does not match the size of this (" << ths.size()
                                         << ")!");
      ScalarType result = 0;
      for (size_t ii = 0; ii < ths.size(); ++ii)
        result += conj(ths.get_unchecked_ref(ii)) * other.get_unchecked_ref(ii);
      return result;
    }
  }; // struct dot_helper<true, ...>

  template <bool anything>
  struct dot_helper<false, anything>
  {
    static ScalarType compute(const derived_type& ths, const derived_type& other)
    {
      if (other.size() != ths.size())
        DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                   "The size of other (" << other.size() << ") does not match the size of this (" << ths.size()
                                         << ")!");
      ScalarType result = 0;
      for (size_t ii = 0; ii < ths.size(); ++ii)
        result += ths.get_unchecked_ref(ii) * other.get_unchecked_ref(ii);
      return result;
    }
  }; // struct dot_helper<false, ...>

public:
  /**
   *  \brief  Computes the scalar products between two vectors.
   *  \param  other The second factor.
   *  \return The scalar product.
   *  \note   If you override this method please use exceptions instead of assertions (for the python bindings).
   */
  virtual ScalarType dot(const derived_type& other) const
  {
    return dot_helper<>::compute(this->as_imp(), other);
  }

  /**
   *  \brief  The l1-norm of the vector.
   *  \return The l1-norm of the vector.
   *  \note   If you override this method please use exceptions instead of assertions (for the python bindings).
   */
  virtual RealType l1_norm() const
  {
    using std::abs;
    RealType result = 0;
    for (size_t ii = 0; ii < size(); ++ii)
      result += abs(get_unchecked_ref(ii));
    return result;
  } // ... l1_norm(...)

  /**
   *  \brief  The l2-norm of the vector.
   *  \return The l2-norm of the vector.
   *  \note   If you override this method please use exceptions instead of assertions (for the python bindings).
   */
  virtual RealType l2_norm() const
  {
    using std::abs;
    using std::sqrt;
    return sqrt(abs(dot(this->as_imp(*this)))); // std::abs is only needed for the right return type:
    // v.dot(v) should always be a ScalarType with zero imaginary part
  }

  /**
   *  \brief  The l-infintiy-norm of the vector.
   *  \return The l-infintiy-norm of the vector.
   *  \note   If you override this method please use exceptions instead of assertions (for the python bindings).
   */
  virtual RealType sup_norm() const
  {
    return amax().second;
  }

  virtual ScalarType standard_deviation() const
  {
    using std::pow;
    using std::sqrt;
    const ScalarType mu = mean();
    ScalarType sigma = 0.0;
    for (const auto& x_i : *this)
      sigma += pow(x_i - mu, 2);
    sigma /= size();
    return sqrt(sigma);
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
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    if (result.size() != size())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of result (" << result.size() << ") does not match the size of this (" << size() << ")!");
    for (size_t ii = 0; ii < size(); ++ii)
      result.set_entry(ii, get_unchecked_ref(ii) + other.get_unchecked_ref(ii));
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
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    for (size_t ii = 0; ii < size(); ++ii)
      add_to_entry(ii, other.get_unchecked_ref(ii));
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
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    if (result.size() != size())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of result (" << result.size() << ") does not match the size of this (" << size() << ")!");
    for (size_t ii = 0; ii < size(); ++ii)
      result.set_entry(ii, get_unchecked_ref(ii) - other.get_unchecked_ref(ii));
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
    static const ScalarType neg_one(-1);
    if (other.size() != size())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "The size of other (" << other.size() << ") does not match the size of this (" << size() << ")!");
    for (size_t ii = 0; ii < size(); ++ii)
      add_to_entry(ii, neg_one * other.get_unchecked_ref(ii));
  } // ... isub(...)

  using BaseType::operator*;

  /**
   *  \brief  Computes the scalar products between this and another vector.
   *  \param  other The second factor.
   *  \return The scalar product.
   *  \see dot()
   */
  virtual ScalarType operator*(const derived_type& other) const
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

  iterator begin()
  {
    this->as_imp().ensure_uniqueness();
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

  operator std::vector<ScalarType>() const
  {
    std::vector<ScalarType> ret(dim());
    for (size_t ii = 0; ii < dim(); ++ii)
      ret[ii] = this->operator[](ii);
    return ret;
  }

private:
  template <class T, class S>
  friend std::ostream& operator<<(std::ostream& /*out*/, const VectorInterface<T, S>& /*vector*/);
}; // class VectorInterface


namespace internal {


template <class VectorImp>
struct VectorAbstractionBase
{
  static const bool is_vector = LA::is_vector<VectorImp>::value;

  static const bool has_static_size = false;

  static const size_t static_size = std::numeric_limits<size_t>::max();

  static const bool is_contiguous = true;

  typedef typename std::conditional<is_vector, VectorImp, void>::type VectorType;
  typedef typename std::conditional<is_vector, typename VectorImp::ScalarType, void>::type ScalarType;
  typedef typename std::conditional<is_vector, typename VectorImp::RealType, void>::type RealType;
  typedef ScalarType S;
  typedef RealType R;

  template <size_t SIZE = static_size, class Field = ScalarType>
  using VectorTypeTemplate = typename std::conditional<is_vector, VectorImp, void>::type;

  template <size_t SIZE = static_size>
  static inline typename std::enable_if<is_vector, VectorType>::type
  create(const size_t sz, const ScalarType& val = Common::suitable_default<ScalarType>::value())
  {
    return VectorType(sz, val);
  }

  static inline ScalarType get_entry(const VectorType& vector, const size_t ii)
  {
    return vector.get_entry(ii);
  }

  static inline void set_entry(VectorType& vector, const size_t ii, const ScalarType& val)
  {
    vector.set_entry(ii, val);
  }

  static inline void add_to_entry(VectorType& vector, const size_t ii, const ScalarType& val)
  {
    vector.add_to_entry(ii, val);
  }

  static inline ScalarType* data(VectorType& vec)
  {
    return &(vec[0]);
  }

  static inline const ScalarType* data(const VectorType& vec)
  {
    return &(vec[0]);
  }
}; // struct VectorAbstractionBase

template <class XtLaVectorImp, size_t size>
struct FieldVectorToLaVector
{
  typedef XtLaVectorImp LaVectorType;
  typedef Dune::FieldVector<typename LaVectorType::ScalarType, size> FieldVectorType;

  static LaVectorType convert(const FieldVectorType& in)
  {
    LaVectorType out(size);
    for (size_t ii = 0; ii < size; ++ii)
      out.set_entry(ii, in[ii]);
    return out;
  }

  static FieldVectorType convert_back(const LaVectorType& in)
  {
    FieldVectorType out;
    for (size_t ii = 0; ii < size; ++ii)
      out[ii] = in.get_entry(ii);
    return out;
  }
};


} // namespace internal


template <class T, class S>
std::ostream& operator<<(std::ostream& out, const VectorInterface<T, S>& vector)
{
  out << "[";
  const size_t sz = vector.size();
  if (sz > 0) {
    out << vector.get_entry(0);
    for (size_t ii = 1; ii < sz; ++ii)
      out << "\n " << vector.get_entry(ii);
  } else
    out << " ";
  out << "]";
  return out;
} // ... operator<<(...)


} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_CONTAINER_VECTOR_INTERFACE_HH
