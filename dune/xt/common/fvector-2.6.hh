// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   René Fritze     (2014 - 2019)
//   Tim Keil        (2018)
//   Tobias Leibner  (2014, 2016 - 2019)

#ifndef DUNE_XT_COMMON_FVECTOR_26_HH
#define DUNE_XT_COMMON_FVECTOR_26_HH

#include <functional>
#include <initializer_list>
#include <numeric>
#include <type_traits>
#include <vector>

#include <boost/functional/hash.hpp>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/promotiontraits.hh>

#include <dune/xt/common/debug.hh>
#include <dune/xt/common/densevector.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/vector.hh>
#include <dune/xt/common/type_traits.hh>

namespace Dune {
namespace XT::Common {


/**
 * \todo We need to implement all operators from the base which return the base, to rather return ourselfes!
 */
template <class K, int SIZE>
class FieldVector : public Dune::FieldVector<K, SIZE>
{
  static_assert(SIZE >= 0, "Really?");

  typedef Dune::FieldVector<K, SIZE> BaseType;
  typedef FieldVector<K, SIZE> ThisType;

public:
  FieldVector(const K& kk = suitable_default<K>::value())
    : BaseType(kk)
  {}

  FieldVector(const ThisType& other) = default;

  FieldVector(const BaseType& other)
    : BaseType(other)
  {}

  /* FieldMatrix< K, 1, 1 > is convertible to K, which in turn is convertible to FieldVector< K, 1 >. Without the
   * following constructor, this leads to an "ambiguous constructor" error (candidates are copy constructor and
   * constructor taking a K) */
  template <class Type = K>
  FieldVector(const typename std::enable_if<SIZE == 1 && std::is_same<K, Type>::value,
                                            typename Dune::FieldMatrix<K, 1, 1>>::type& mat)
    : BaseType(mat[0][0])
  {}

  FieldVector(const std::vector<K>& vec)
    : BaseType()
  {
#ifndef NDEBUG
    if (vec.size() != SIZE)
      DUNE_THROW(Exceptions::wrong_input_given,
                 "You are trying to construct a FieldVector< ..., " << SIZE << " > (of "
                                                                    << "static size) from a vector of size "
                                                                    << vec.size() << "!");
#endif // NDEBUG
    for (size_t ii = 0; ii < SIZE; ++ii)
      this->operator[](ii) = vec[ii];
  } // FieldVector(...)

  FieldVector(std::initializer_list<K> list)
    : BaseType()
  {
#ifndef NDEBUG
    if (list.size() != SIZE)
      DUNE_THROW(Exceptions::wrong_input_given,
                 "You are trying to construct a FieldVector< ..., " << SIZE << " > (of "
                                                                    << "static size) from a list of size "
                                                                    << list.size() << "!");
#endif // NDEBUG
    size_t ii = 0;
    for (auto element : list)
      this->operator[](ii++) = element;
  } // FieldVector(...)

  template <class C>
  FieldVector(
      const DenseVector<C>& x,
      typename std::enable_if<IsFieldVectorSizeCorrect<C, SIZE>::value
                              && std::is_convertible<typename DenseVector<C>::value_type, K>::value>::type* /*dummy*/
      = nullptr)
    : BaseType(x)
  {}

  operator std::vector<K>() const
  {
    std::vector<K> ret(SIZE);
    for (size_t ii = 0; ii < SIZE; ++ii)
      ret[ii] = this->operator[](ii);
    return ret;
  }

  operator Dune::FieldMatrix<K, 1, SIZE>() const
  {
    Dune::FieldMatrix<K, 1, SIZE> ret;
    ret[0] = *this;
    return ret;
  }

  template <int S = SIZE>
  operator typename std::enable_if<(S == SIZE) && (SIZE != 1), Dune::FieldMatrix<K, S, 1>>::type() const
  {
    Dune::FieldMatrix<K, SIZE, 1> ret;
    for (size_t ii = 0; ii < SIZE; ++ii)
      ret[ii][0] = this->operator[](ii);
    return ret;
  }

  operator std::array<K, SIZE>() const
  {
    std::array<K, SIZE> ret;
    for (size_t ii = 0; ii < SIZE; ++ii)
      ret[ii] = this->operator[](ii);
    return ret;
  }

  template <class K_Other>
  typename multiplication_promotion<K, K_Other>::type operator*(const Dune::FieldVector<K_Other, SIZE>& other) const
  {
    typename multiplication_promotion<K, K_Other>::type result(0.);
    for (size_t ii = 0; ii < SIZE; ++ii)
      result += (*this)[ii] * other[ii];
    return result;
  }

  template <int R>
  typename std::enable_if<SIZE != 1 && R == SIZE, K>::type operator*(const Dune::FieldMatrix<K, R, 1>& mat) const
  {
    K ret = 0;
    for (size_t ii = 0; ii < SIZE; ++ii)
      ret += (*this)[ii] * mat[ii][0];
    return ret;
  }

  //! This op is not redundant
  ThisType operator*(const K& scalar) const
  {
    ThisType ret(*this);
    ret *= scalar;
    return ret;
  }

  ThisType operator/(const K& scalar) const
  {
    ThisType ret(*this);
    ret /= scalar;
    return ret;
  }
}; // class FieldVector


template <class K, size_t block_num, size_t size_block>
class BlockedFieldVector
{
  using ThisType = BlockedFieldVector;

public:
  static constexpr size_t num_blocks = block_num;
  static constexpr size_t block_size = size_block;
  static constexpr size_t static_size = num_blocks * block_size;
  using VectorType = Dune::FieldVector<K, static_size>;
  using XtVectorType = FieldVector<K, static_size>;
  using BlockType = FieldVector<K, block_size>;

  BlockedFieldVector(const K& val = K(0.))
    : backend_(BlockType(val))
  {}

  template <class OtherVectorType>
  BlockedFieldVector(const OtherVectorType& other,
                     std::enable_if_t<is_vector<OtherVectorType>::value, int> /*dummy*/ = 0)
  {
    *this = other;
  }

  BlockedFieldVector(const BlockType& block)
    : backend_(block)
  {}

  BlockedFieldVector(const std::initializer_list<K>& init_list)
  {
    assert(init_list.size() == static_size || init_list.size() == num_blocks);
    if (init_list.size() == static_size) {
      *this = VectorType(init_list);
    } else {
      size_t jj = 0;
      for (const auto& block_init : init_list)
        block(jj++) = block_init;
    }
  }

  template <class OtherVectorType>
  std::enable_if_t<is_vector<OtherVectorType>::value, ThisType>& operator=(const OtherVectorType& other)
  {
    assert(other.size() == static_size);
    for (size_t jj = 0; jj < num_blocks; ++jj)
      for (size_t ii = 0; ii < block_size; ++ii)
        backend_[jj][ii] = other[jj * block_size + ii];
    return *this;
  }

  size_t size() const
  {
    return static_size;
  }

  K& operator[](const size_t ii)
  {
    assert(ii < static_size);
    return backend_[ii / block_size][ii % block_size];
  }

  const K& operator[](const size_t ii) const
  {
    assert(ii < static_size);
    return backend_[ii / block_size][ii % block_size];
  }

  K& get_entry(const size_t jj, const size_t ii)
  {
    assert(jj < num_blocks);
    assert(ii < block_size);
    return backend_[jj][ii];
  }

  const K& get_entry(const size_t jj, const size_t ii) const
  {
    assert(jj < num_blocks);
    assert(ii < block_size);
    return backend_[jj][ii];
  }

  BlockType& block(const size_t jj)
  {
    assert(jj < num_blocks);
    return backend_[jj];
  }

  const BlockType& block(const size_t jj) const
  {
    assert(jj < num_blocks);
    return backend_[jj];
  }

  K one_norm() const
  {
    K ret(0);
    for (const auto& block_vec : backend_)
      for (const auto& entry : block_vec)
        ret += std::abs(entry);
    return ret;
  }

  K two_norm() const
  {
    return std::sqrt(two_norm2());
  }

  K two_norm2() const
  {
    K ret(0);
    for (const auto& block_vec : backend_)
      for (const auto& entry : block_vec)
        ret += std::pow(entry, 2);
    return ret;
  }

  operator XtVectorType() const
  {
    XtVectorType ret(0.);
    for (size_t jj = 0; jj < num_blocks; ++jj)
      for (size_t ii = 0; ii < block_size; ++ii)
        ret[jj * block_size + ii] = backend_[jj][ii];
    return ret;
  }

  K* data()
  {
    return &(backend_[0][0]);
  }

  const K* data() const
  {
    return &(backend_[0][0]);
  }

  K* begin()
  {
    return data();
  }

  const K* begin() const
  {
    return data();
  }

  K* end()
  {
    return &(backend_[num_blocks - 1][block_size]);
  }

  const K* end() const
  {
    return &(backend_[num_blocks - 1][block_size]);
  }

  ThisType& operator*=(const K& val)
  {
    for (size_t jj = 0; jj < num_blocks; ++jj)
      backend_[jj] *= val;
    return *this;
  }

  ThisType operator*(const K& val) const
  {
    auto ret = *this;
    ret *= val;
    return ret;
  }

  K operator*(const ThisType& other) const
  {
    return std::inner_product(begin(), end(), other.begin(), 0.);
  }

  ThisType& operator+=(const ThisType& other)
  {
    backend_ += other.backend_;
    return *this;
  }

  ThisType operator+(const ThisType& other) const
  {
    auto ret = *this;
    ret += other;
    return ret;
  }

  ThisType& operator-=(const ThisType& other)
  {
    backend_ -= other.backend_;
    return *this;
  }

  ThisType operator-(const ThisType& other) const
  {
    auto ret = *this;
    ret -= other;
    return ret;
  }

private:
  FieldVector<BlockType, num_blocks> backend_;
};


//! this allows to set the init value of the FieldVector at compile time
template <class K, int SIZE, K value>
class ValueInitFieldVector : public Dune::XT::Common::FieldVector<K, SIZE>
{
  typedef Dune::XT::Common::FieldVector<K, SIZE> BaseType;

public:
  ValueInitFieldVector()
    : BaseType(value)
  {}
}; // class ValueInitFieldVector


//! struct to be used as comparison function e.g. in a std::map<FieldVector<...>, ..., FieldVectorLess>
struct FieldVectorLess
{
  template <class FieldType, int dimDomain>
  bool operator()(const Dune::FieldVector<FieldType, dimDomain>& a,
                  const Dune::FieldVector<FieldType, dimDomain>& b) const
  {
    for (size_t dd = 0; dd < dimDomain; ++dd) {
      if (a[dd] < b[dd])
        return true;
      else if (a[dd] > b[dd])
        return false;
    }
    return false;
  }
};

struct FieldVectorFloatLess
{
  template <class FieldType, int dimDomain>
  bool operator()(const Dune::FieldVector<FieldType, dimDomain>& a,
                  const Dune::FieldVector<FieldType, dimDomain>& b) const
  {
    for (size_t dd = 0; dd < dimDomain; ++dd) {
      if (XT::Common::FloatCmp::lt(a[dd], b[dd]))
        return true;
      else if (XT::Common::FloatCmp::gt(a[dd], b[dd]))
        return false;
    }
    return false;
  }
};


//! Specialization of VectorAbstraction for Dune::XT::Common::FieldVector
template <class K, int SIZE>
struct VectorAbstraction<Dune::XT::Common::FieldVector<K, SIZE>>
  : public internal::VectorAbstractionBase<Dune::XT::Common::FieldVector<K, SIZE>, K>
  , public internal::HasSubscriptOperatorForVectorAbstraction<Dune::XT::Common::FieldVector<K, SIZE>,
                                                              typename Dune::FieldTraits<K>::field_type>
{
  static const bool has_static_size = true;
  static const size_t static_size = SIZE;
  static const bool is_contiguous = true;

  template <size_t SZ = SIZE, class Field = K>
  using VectorTypeTemplate = Dune::XT::Common::FieldVector<Field, SZ>;

  template <size_t SZ = SIZE>
  static inline VectorTypeTemplate<SZ> create(const size_t sz, const K& val = suitable_default<K>::value())
  {
    if (sz != SZ)
      DUNE_THROW(Exceptions::wrong_input_given,
                 "You are trying to construct a FieldVector< ..., " << SZ << " > (of "
                                                                    << "static size) with " << sz << " elements!");
    return VectorTypeTemplate<SZ>(val);
  }
};


//! Specialization of VectorAbstraction for Dune::XT::Common::BlockedFieldVector
template <class K, size_t num_blocks, size_t block_size>
struct VectorAbstraction<Dune::XT::Common::BlockedFieldVector<K, num_blocks, block_size>>
  : public internal::VectorAbstractionBase<Dune::XT::Common::BlockedFieldVector<K, num_blocks, block_size>, K>
  , public internal::HasSubscriptOperatorForVectorAbstraction<
        Dune::XT::Common::BlockedFieldVector<K, num_blocks, block_size>,
        typename Dune::FieldTraits<K>::field_type>
{
  using VectorType = Dune::XT::Common::BlockedFieldVector<K, num_blocks, block_size>;
  static constexpr bool has_static_size = true;
  static constexpr size_t static_size = VectorType::static_size;
  static constexpr bool is_contiguous = true;

  template <size_t SZ = static_size>
  static inline VectorType create(const size_t sz, const K& val = suitable_default<K>::value())
  {
    static_assert(SZ == static_size, "Creation of Vector with different size not implemented!");
    if (sz != SZ)
      DUNE_THROW(Exceptions::wrong_input_given,
                 "You are trying to construct a FieldVector< ..., " << SZ << " > (of "
                                                                    << "static size) with " << sz << " elements!");
    return VectorType(val);
  }
};


template <class V>
typename std::enable_if<
    is_vector<V>::value && VectorAbstraction<V>::has_static_size,
    std::unique_ptr<FieldVector<typename VectorAbstraction<V>::S, VectorAbstraction<V>::static_size>>>::type
make_field_container_ptr(const V& vec)
{
  auto ret = std::make_unique<FieldVector<typename VectorAbstraction<V>::S, VectorAbstraction<V>::static_size>>;
  for (size_t ii = 0; ii < ret->size(); ++ii)
    (*ret)[ii] = vec[ii];
  return std::move(ret);
}


template <class V>
typename std::enable_if<is_vector<V>::value && VectorAbstraction<V>::has_static_size,
                        FieldVector<typename VectorAbstraction<V>::S, VectorAbstraction<V>::static_size>>::type
make_field_container(const V& vec)
{
  FieldVector<typename VectorAbstraction<V>::S, VectorAbstraction<V>::static_size> ret;
  for (size_t ii = 0; ii < ret.size(); ++ii)
    ret[ii] = vec[ii];
  return ret;
}


template <class K, int SIZE>
FieldVector<K, SIZE> make_field_container(Dune::FieldVector<K, SIZE>&& vec)
{
  return std::move(vec);
}


template <class K, int SIZE>
FieldVector<K, SIZE> make_field_container(FieldVector<K, SIZE>&& vec)
{
  return std::move(vec);
}

template <class K>
FieldVector<K, 3> cross_product(const FieldVector<K, 3>& u, const FieldVector<K, 3>& v)
{
  FieldVector<K, 3> ret;
  ret[0] = u[1] * v[2] - u[2] * v[1];
  ret[1] = u[2] * v[0] - u[0] * v[2];
  ret[2] = u[0] * v[1] - u[1] * v[0];
  return ret;
}


namespace internal {


/**
 * Most likely, you do not want to use this class directly, but \sa hstack instead.
 */
template <class K>
struct hstack_decay
{
  static_assert(is_arithmetic<K>::value);
  using type = K;
};

template <class K, int SIZE>
struct hstack_decay<Dune::FieldVector<K, SIZE>>
{
  using type = Dune::XT::Common::FieldVector<K, SIZE>;
};

template <class K, int SIZE>
struct hstack_decay<const Dune::FieldVector<K, SIZE>&> : public hstack_decay<Dune::FieldVector<K, SIZE>>
{};

template <class K, int SIZE>
struct hstack_decay<Dune::XT::Common::FieldVector<K, SIZE>> : public hstack_decay<Dune::FieldVector<K, SIZE>>
{};

template <class K, int SIZE>
struct hstack_decay<const Dune::XT::Common::FieldVector<K, SIZE>&> : public hstack_decay<Dune::FieldVector<K, SIZE>>
{};


/**
 * Most likely, you do not want to use this, but \sa hstack instead.
 */
template <class T>
using hstack_decay_t = typename hstack_decay<T>::type;


/**
 * Most likely, you do not want to use this class directly, but \sa hstack instead.
 */
template <class... Vectors>
struct hstack_helper;

template <class L, class R, class... Vectors>
struct hstack_helper<L, R, Vectors...>
{
  using type =
      typename hstack_helper<typename hstack_helper<hstack_decay_t<L>, hstack_decay_t<R>>::type, Vectors...>::type;
};

template <class K, int SIZE>
struct hstack_helper<Dune::FieldVector<K, SIZE>>
{
  using type = hstack_decay_t<Dune::FieldVector<K, SIZE>>;
};

template <class K, int SIZE>
struct hstack_helper<Dune::XT::Common::FieldVector<K, SIZE>> : public hstack_helper<Dune::FieldVector<K, SIZE>>
{};

template <class KL, class KR, int r>
struct hstack_helper<KL, Dune::FieldVector<KR, r>>
{
  using type = hstack_decay_t<Dune::FieldVector<typename PromotionTraits<KL, KR>::PromotedType, 1 + r>>;
};

template <class KL, class KR, int r>
struct hstack_helper<KL, Dune::XT::Common::FieldVector<KR, r>> : public hstack_helper<KL, Dune::FieldVector<KR, r>>
{};

template <class KL, int l, class KR>
struct hstack_helper<Dune::FieldVector<KL, l>, KR>
{
  using type = hstack_decay_t<Dune::FieldVector<typename PromotionTraits<KL, KR>::PromotedType, l + 1>>;
};

template <class KL, int l, class KR>
struct hstack_helper<Dune::XT::Common::FieldVector<KL, l>, KR> : public hstack_helper<Dune::FieldVector<KL, l>, KR>
{};

template <class KL, int l, class KR, int r>
struct hstack_helper<Dune::FieldVector<KL, l>, Dune::FieldVector<KR, r>>
{
  using type = hstack_decay_t<Dune::FieldVector<typename PromotionTraits<KL, KR>::PromotedType, l + r>>;
};

template <class KL, int l, class KR, int r>
struct hstack_helper<Dune::XT::Common::FieldVector<KL, l>, Dune::FieldVector<KR, r>>
  : public hstack_helper<Dune::FieldVector<KL, l>, Dune::FieldVector<KR, r>>
{};

template <class KL, int l, class KR, int r>
struct hstack_helper<Dune::FieldVector<KL, l>, Dune::XT::Common::FieldVector<KR, r>>
  : public hstack_helper<Dune::FieldVector<KL, l>, Dune::FieldVector<KR, r>>
{};

template <class KL, int l, class KR, int r>
struct hstack_helper<Dune::XT::Common::FieldVector<KL, l>, Dune::XT::Common::FieldVector<KR, r>>
  : public hstack_helper<Dune::FieldVector<KL, l>, Dune::FieldVector<KR, r>>
{};

// not sure why this is required, would've hoped the decay above would take care of it
template <class KL, int l, class KR, int r>
struct hstack_helper<Dune::FieldVector<KL, l>, const Dune::FieldVector<KR, r>&>
  : public hstack_helper<Dune::FieldVector<KL, l>, Dune::FieldVector<KR, r>>
{};

// not sure why this is required, would've hoped the decay above would take care of it
template <class KL, int l, class KR, int r>
struct hstack_helper<Dune::FieldVector<KL, l>, const Dune::XT::Common::FieldVector<KR, r>&>
  : public hstack_helper<Dune::FieldVector<KL, l>, Dune::FieldVector<KR, r>>
{};


} // namespace internal


/// \brief Specialization of \sa hstack to stop the recursion.
template <class K, int SIZE>
Dune::XT::Common::FieldVector<K, SIZE> hstack(Dune::FieldVector<K, SIZE>&& vec)
{
  return std::move(vec);
}

/// \brief Specialization of \sa hstack to stop the recursion.
template <class K, int SIZE>
Dune::XT::Common::FieldVector<K, SIZE> hstack(Dune::XT::Common::FieldVector<K, SIZE>&& vec)
{
  return std::move(vec);
}

/// \brief Specialization of \sa hstack to stop the recursion.
template <class K, int SIZE>
Dune::XT::Common::FieldVector<K, SIZE> hstack(const Dune::FieldVector<K, SIZE>& vec)
{
  return vec;
}

/// \brief Specialization of \sa hstack to stop the recursion.
template <class K, int SIZE>
Dune::XT::Common::FieldVector<K, SIZE> hstack(const Dune::XT::Common::FieldVector<K, SIZE>& vec)
{
  return vec;
}


/// \brief Specialization of \sa hstack for two arguments.
template <class KL, class KR, int r>
typename std::enable_if<
    is_arithmetic<KL>::value && !is_field_vector<KL>::value,
    typename internal::hstack_helper<internal::hstack_decay_t<KL>,
                                     internal::hstack_decay_t<Dune::FieldVector<KR, r>>>::type>::type
hstack(const KL& left, const Dune::FieldVector<KR, r>& right)
{
  typename internal::hstack_helper<internal::hstack_decay_t<KL>,
                                   internal::hstack_decay_t<Dune::FieldVector<KR, r>>>::type ret;
  ret[0] = left;
  for (size_t ii = 0; ii < r; ++ii)
    ret[ii + 1] = right[ii];
  return ret;
} // ... hstack<1, r>(...)

/// \brief Specialization of \sa hstack for two arguments.
template <class KL, class KR, int r>
typename std::enable_if<
    is_arithmetic<KL>::value && !is_field_vector<KL>::value,
    typename internal::hstack_helper<internal::hstack_decay_t<KL>,
                                     internal::hstack_decay_t<Dune::FieldVector<KR, r>>>::type>::type
hstack(const KL& left, const Dune::XT::Common::FieldVector<KR, r>& right)
{
  return hstack(left, static_cast<const Dune::FieldVector<KR, r>&>(right));
}

/// \brief Specialization of \sa hstack for two arguments.
template <class KL, int l, class KR>
typename std::enable_if<is_arithmetic<KR>::value && !is_field_vector<KR>::value,
                        typename internal::hstack_helper<internal::hstack_decay_t<Dune::FieldVector<KL, l>>,
                                                         internal::hstack_decay_t<KR>>::type>::type
hstack(const Dune::FieldVector<KL, l>& left, const KR& right)
{
  typename internal::hstack_helper<internal::hstack_decay_t<Dune::FieldVector<KL, l>>,
                                   internal::hstack_decay_t<KR>>::type ret;
  for (size_t ii = 0; ii < l; ++ii)
    ret[ii] = left[ii];
  ret[l] = right;
  return ret;
} // ... hstack<l, 1>(...)

/// \brief Specialization of \sa hstack for two arguments.
template <class KL, int l, class KR>
typename std::enable_if<is_arithmetic<KR>::value && !is_field_vector<KR>::value,
                        typename internal::hstack_helper<internal::hstack_decay_t<Dune::FieldVector<KL, l>>,
                                                         internal::hstack_decay_t<KR>>::type>::type
hstack(const Dune::XT::Common::FieldVector<KL, l>& left, const KR& right)
{
  return hstack(static_cast<const Dune::FieldVector<KL, l>&>(left), right);
}

/// \brief Specialization of \sa hstack for two arguments.
template <class KL, int l, class KR, int r>
typename internal::hstack_helper<internal::hstack_decay_t<Dune::FieldVector<KL, l>>,
                                 internal::hstack_decay_t<Dune::FieldVector<KR, r>>>::type
hstack(const Dune::FieldVector<KL, l>& left, const Dune::FieldVector<KR, r>& right)
{
  typename internal::hstack_helper<internal::hstack_decay_t<Dune::FieldVector<KL, l>>,
                                   internal::hstack_decay_t<Dune::FieldVector<KR, r>>>::type ret;
  size_t ii = 0;
  for (size_t jj = 0; jj < l; ++jj, ++ii)
    ret[ii] = left[jj];
  for (size_t jj = 0; jj < r; ++jj, ++ii)
    ret[ii] = right[jj];
  return ret;
} // ... hstack<l, r>(...)

/// \brief Specialization of \sa hstack for two arguments.
template <class KL, int l, class KR, int r>
typename internal::hstack_helper<internal::hstack_decay_t<Dune::FieldVector<KL, l>>,
                                 internal::hstack_decay_t<Dune::FieldVector<KR, r>>>::type
hstack(const Dune::XT::Common::FieldVector<KL, l>& left, const Dune::FieldVector<KR, r>& right)
{
  return hstack(static_cast<const Dune::FieldVector<KL, l>&>(left), right);
}

/// \brief Specialization of \sa hstack for two arguments.
template <class KL, int l, class KR, int r>
typename internal::hstack_helper<internal::hstack_decay_t<Dune::FieldVector<KL, l>>,
                                 internal::hstack_decay_t<Dune::FieldVector<KR, r>>>::type
hstack(const Dune::FieldVector<KL, l>& left, const Dune::XT::Common::FieldVector<KR, r>& right)
{
  return hstack(left, static_cast<const Dune::FieldVector<KR, r>&>(right));
}

/// \brief Specialization of \sa hstack for two arguments.
template <class KL, int l, class KR, int r>
typename internal::hstack_helper<internal::hstack_decay_t<Dune::FieldVector<KL, l>>,
                                 internal::hstack_decay_t<Dune::FieldVector<KR, r>>>::type
hstack(const Dune::XT::Common::FieldVector<KL, l>& left, const Dune::XT::Common::FieldVector<KR, r>& right)
{
  return hstack(static_cast<const Dune::FieldVector<KL, l>&>(left),
                static_cast<const Dune::FieldVector<KR, r>&>(right));
}


/**
 * \brief Stacks an arbitrary number of numbers and FieldVectors horizontally into one appropriate FieldVector.
 * \note  At least one argument must be a FieldVector!
 */
template <class L, class R, class... Vectors>
typename std::enable_if<
    is_field_vector<L>::value || is_field_vector<R>::value,
    typename internal::hstack_helper<internal::hstack_decay_t<L>, internal::hstack_decay_t<R>, Vectors...>::type>::type
hstack(const L& left, const R& right, Vectors&&... vectors)
{
  return hstack(hstack(left, right), std::forward<Vectors>(vectors)...);
}


} // namespace XT::Common


template <class K, int SIZE>
typename std::enable_if<SIZE != 1, K>::type operator*(const Dune::FieldVector<K, SIZE>& vec,
                                                      const Dune::FieldMatrix<K, SIZE, 1>& mat)
{
  K ret = 0;
  for (size_t ii = 0; ii < SIZE; ++ii)
    ret += vec[ii] * mat[ii][0];
  return ret;
}

template <class K, int SIZE>
struct FieldTraits<XT::Common::FieldVector<K, SIZE>>
{
  using field_type = typename FieldTraits<K>::field_type;
  using real_type = typename FieldTraits<K>::real_type;
};


} // namespace Dune
namespace std {


/**
 * \sa http://en.cppreference.com/w/cpp/utility/hash
 * \sa http://www.boost.org/doc/libs/1_55_0/doc/html/hash/combine.html
 */
template <class K, int SIZE>
struct hash<Dune::FieldVector<K, SIZE>>
{
  typedef Dune::FieldVector<K, SIZE> argument_type;
  using result_type = std::size_t;

  result_type operator()(argument_type const& s) const noexcept
  {
    std::size_t seed = 0;
    for (size_t ii = 0; ii < s.size(); ++ii)
      boost::hash_combine(seed, s[ii]);
    return seed;
  }
}; // struct hash<Dune::FieldVector<...>>


/**
 * \sa http://en.cppreference.com/w/cpp/utility/hash
 * \sa http://www.boost.org/doc/libs/1_55_0/doc/html/hash/combine.html
 */
template <class K, int SIZE>
struct hash<Dune::XT::Common::FieldVector<K, SIZE>>
{
  typedef Dune::XT::Common::FieldVector<K, SIZE> argument_type;
  using result_type = std::size_t;

  result_type operator()(argument_type const& s) const noexcept
  {
    std::size_t seed = 0;
    for (size_t ii = 0; ii < s.size(); ++ii)
      boost::hash_combine(seed, s[ii]);
    return seed;
  }
}; // struct hash<Dune::XT::Common::FieldVector<...>>


} // namespace std

#endif // DUNE_XT_COMMON_FVECTOR_26_HH
