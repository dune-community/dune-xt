// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Barbara Verfürth (2015)
//   Felix Schindler  (2012 - 2017, 2019)
//   René Fritze      (2012 - 2019)
//   Sven Kaulmann    (2013)
//   Tobias Leibner   (2014 - 2020)

#ifndef DUNE_XT_COMMON_VECTOR_HH
#define DUNE_XT_COMMON_VECTOR_HH

#include <vector>
#include <ostream>

#include <dune/common/densevector.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/fvector.hh>
#include <dune/common/ftraits.hh>
#include <dune/common/streamoperators.hh> // for operator<<(stream, std::array<>)

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/numeric_cast.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/common/debug.hh>

namespace Dune {
namespace XT {
namespace Common {


namespace internal {


/**
 * \brief A base class to derive from when specializing \sa VectorAbstraction.
 */
template <class Vector, class K>
struct VectorAbstractionBase
{
  typedef Vector VectorType;
  using ScalarType = typename Dune::FieldTraits<K>::field_type;
  using RealType = typename Dune::FieldTraits<K>::real_type;
  using S = ScalarType;
  using R = RealType;

  static const constexpr bool is_vector = true;
  static const constexpr bool has_static_size = false;
  static const constexpr size_t static_size = std::numeric_limits<size_t>::max();
  static const constexpr bool has_ostream = false;
  static const constexpr bool is_contiguous = false;

  template <size_t SIZE = static_size, class FieldType = ScalarType>
  using VectorTypeTemplate = Vector;
}; // struct VectorAbstractionBase


/**
 * \brief A base class to derive from when specializing \sa VectorAbstraction.
 */
template <class VectorType, class ScalarType>
struct HasSubscriptOperatorForVectorAbstraction
{
  static inline ScalarType get_entry(const VectorType& vector, const size_t ii)
  {
    DXT_ASSERT(ii < vector.size());
    return vector[ii];
  }

  static inline void set_entry(VectorType& vector, const size_t ii, const ScalarType& val)
  {
    DXT_ASSERT(ii < vector.size());
    vector[ii] = val;
  }

  static inline void add_to_entry(VectorType& vector, const size_t ii, const ScalarType& val)
  {
    DXT_ASSERT(ii < vector.size());
    vector[ii] += val;
  }

  static inline ScalarType* data(VectorType& vec)
  {
    return &(vec[0]);
  }

  static inline const ScalarType* data(const VectorType& vec)
  {
    return &(vec[0]);
  }
};


} // namespace internal


/**
 * \brief Traits to uniformly handle dense (mathematical) vectors.
 *
 *        If you want your vector class to benefit from the functionality in this header you have to manually
 *        specify a specialization of this class in your code with is_vector defined to true and the appropriate
 *        static methods implemented and members defined (see the specializations below).
 * \sa internal::VectorAbstractionBase
 * \sa internal::HasSubscriptOperatorForVectorAbstraction
 */
template <class VecType>
struct VectorAbstraction
{
  using V = std::conditional_t<std::is_same<VecType, void>::value, int, VecType>; // avoid reference to void
  typedef V VectorType;
  typedef typename Dune::FieldTraits<V>::field_type ScalarType;
  typedef typename Dune::FieldTraits<V>::real_type RealType;
  typedef ScalarType S;
  typedef RealType R;

  static const constexpr bool is_vector = false;
  static const constexpr bool has_static_size = false;
  static const constexpr size_t static_size = std::numeric_limits<size_t>::max();
  static const constexpr bool has_ostream = true;
  static const constexpr bool is_contiguous = false;

  template <size_t SIZE = static_size, class FieldType = ScalarType>
  using VectorTypeTemplate = VecType;

  static inline /*VectorType*/ void create(const size_t /*sz*/)
  {
    static_assert(AlwaysFalse<VecType>::value, "Do not call me if is_vector is false!");
  }

  static inline /*VectorType*/ void create(const size_t /*sz*/, const ScalarType& /*val*/)
  {
    static_assert(AlwaysFalse<VecType>::value, "Do not call me if is_vector is false!");
  }

  static inline /*ScalarType*/ void get_entry(const VectorType& /*vec*/, const size_t /*ii*/)
  {
    static_assert(AlwaysFalse<VecType>::value, "Do not call me if is_vector is false!");
  }

  static inline void set_entry(VectorType& /*vec*/, size_t /*ii*/, ScalarType& /*val*/)
  {
    static_assert(AlwaysFalse<VecType>::value, "Do not call me if is_vector is false!");
  }

  static inline ScalarType* data(std::remove_const_t<VectorType>& /*vec*/)
  {
    static_assert(AlwaysFalse<VecType>::value, "Do not call me if is_contiguous is false!");
    return nullptr;
  }

  // using std::enable_if to avoid 'multiple overloads with same signature' for some strange VectorTypes
  template <bool is_vec = is_vector>
  static inline const std::enable_if_t<is_vec, ScalarType>* data(std::add_const_t<VectorType>& /*vec*/)
  {
    static_assert(AlwaysFalse<VecType>::value, "Do not call me if is_contiguous is false!");
    return nullptr;
  }
};


template <class T, class Allocator>
struct VectorAbstraction<std::vector<T, Allocator>>
  : public internal::VectorAbstractionBase<std::vector<T, Allocator>, T>
  , public internal::HasSubscriptOperatorForVectorAbstraction<std::vector<T, Allocator>, T>
{
  static const constexpr bool is_contiguous = true;
  static const constexpr bool static_size = std::numeric_limits<size_t>::max();

  template <size_t SIZE = static_size>
  static inline std::vector<T, Allocator> create(const size_t sz, const T& val = T(0))
  {
    return std::vector<T, Allocator>(sz, val);
  }

  static inline T* data(std::vector<T, Allocator>& vec)
  {
    return vec.data();
  }

  static inline const T* data(const std::vector<T, Allocator>& vec)
  {
    return vec.data();
  }
};

template <class K, size_t SIZE>
struct VectorAbstraction<std::array<K, SIZE>>
  : public internal::VectorAbstractionBase<std::array<K, SIZE>, K>
  , public internal::HasSubscriptOperatorForVectorAbstraction<std::array<K, SIZE>,
                                                              typename Dune::FieldTraits<K>::field_type>
{
  static const constexpr bool has_static_size = true;
  static const constexpr size_t static_size = SIZE;
  static const constexpr bool has_ostream = true; // provided in dune/common/streamoperators.hh
  static const constexpr bool is_contiguous = true;

  template <size_t SZ = SIZE, class FieldType = K>
  using VectorTypeTemplate = std::array<FieldType, SZ>;

  template <size_t SZ = SIZE>
  static inline VectorTypeTemplate<SZ> create(const size_t sz, const K& val = K(0))
  {
    if (sz != SZ)
      DUNE_THROW(Exceptions::shapes_do_not_match, "requested size: " << sz << "\n   actual size: " << int(SZ));
    VectorTypeTemplate<SZ> ret;
    ret.fill(val);
    return ret;
  }

  static inline K* data(std::array<K, SIZE>& vec)
  {
    return vec.data();
  }

  static inline const K* data(const std::array<K, SIZE>& vec)
  {
    return vec.data();
  }
};

template <class K>
struct VectorAbstraction<Dune::DynamicVector<K>>
  : public internal::VectorAbstractionBase<Dune::DynamicVector<K>, K>
  , public internal::HasSubscriptOperatorForVectorAbstraction<Dune::DynamicVector<K>,
                                                              typename Dune::FieldTraits<K>::field_type>
{
  static const constexpr bool is_contiguous = true;
  static const constexpr bool static_size = internal::VectorAbstractionBase<Dune::DynamicVector<K>, K>::static_size;

  template <size_t SIZE = static_size>
  static inline Dune::DynamicVector<K> create(const size_t sz, const K& val = K(0))
  {
    return Dune::DynamicVector<K>(sz, val);
  }
};

template <class K, int SIZE>
struct VectorAbstraction<Dune::FieldVector<K, SIZE>>
  : public internal::VectorAbstractionBase<Dune::FieldVector<K, SIZE>, K>
  , public internal::HasSubscriptOperatorForVectorAbstraction<Dune::FieldVector<K, SIZE>,
                                                              typename Dune::FieldTraits<K>::field_type>
{
  static const constexpr bool has_static_size = true;
  static const constexpr size_t static_size = SIZE;
  static const constexpr bool is_contiguous = true;

  template <size_t SZ = SIZE, class Field = K>
  using VectorTypeTemplate = Dune::FieldVector<Field, SZ>;

  template <size_t SZ = SIZE>
  static inline VectorTypeTemplate<SZ> create(const size_t sz, const K& val = K(0))
  {
    if (sz != SZ)
      DUNE_THROW(Exceptions::shapes_do_not_match, "requested size: " << sz << "\n   actual size: " << int(SZ));
    return VectorTypeTemplate<SZ>(val);
  }
};

template <class V>
struct VectorAbstraction<Dune::DenseVector<V>> : public VectorAbstraction<typename Dune::DenseVector<V>::derived_type>
{};


template <class VectorType, size_t SIZE = VectorAbstraction<VectorType>::static_size>
typename std::enable_if<is_vector<VectorType>::value,
                        typename VectorAbstraction<VectorType>::template VectorTypeTemplate<SIZE>>::type
create(const size_t sz, const typename VectorAbstraction<VectorType>::ScalarType& val = 0)
{
  return VectorAbstraction<VectorType>::template create<SIZE>(sz, val);
}


template <class TargetVectorType, class SourceVectorType>
typename std::enable_if<is_vector<TargetVectorType>::value && is_vector<SourceVectorType>::value,
                        TargetVectorType>::type
zeros_like(const SourceVectorType& source)
{
  return create<TargetVectorType>(source.size(), 0);
}


template <class VectorType>
typename std::enable_if<is_vector<VectorType>::value, VectorType>::type zeros_like(const VectorType& source)
{
  return zeros_like<VectorType, VectorType>(source);
}


template <class VectorType>
std::enable_if_t<is_vector<VectorType>::value, typename VectorAbstraction<VectorType>::R>
sup_norm(const VectorType& vector)
{
  typename VectorAbstraction<VectorType>::R result(0);
  for (size_t ii = 0; ii < vector.size(); ++ii)
    result = std::max(result, std::abs(vector[ii]));
  return result;
}

template <class F>
typename FieldTraits<F>::real_type sup_norm(const DynamicVector<F>& vector)
{
  return vector.infinity_norm();
}

template <class F, int SIZE>
typename FieldTraits<F>::real_type sup_norm(const FieldVector<F, SIZE>& vector)
{
  return vector.infinity_norm();
}


template <class VectorType>
typename std::enable_if<is_vector<VectorType>::value, typename VectorAbstraction<VectorType>::ScalarType*>::type
data(VectorType& source)
{
  return VectorAbstraction<VectorType>::data(source);
}


template <class VectorType>
typename std::enable_if<is_vector<VectorType>::value, const typename VectorAbstraction<VectorType>::ScalarType*>::type
data(const VectorType& source)
{
  return VectorAbstraction<VectorType>::data(source);
}


template <class V, StorageLayout storage_layout = StorageLayout::dense_row_major>
typename std::enable_if<is_vector<V>::value, std::unique_ptr<typename VectorAbstraction<V>::S[]>>::type
serialize(const V& vec)
{
  typedef typename VectorAbstraction<V>::S S;
  static_assert(storage_layout == StorageLayout::dense_row_major || storage_layout == StorageLayout::dense_column_major,
                "You have to select either row_major or column_major as StorageLayout!");
  using Vec = VectorAbstraction<V>;
  const size_t size = Vec::size(vec);
  auto data = std::make_unique<S[]>(size);
  for (size_t ii = 0; ii < size; ++ii)
#ifdef DXT_DISABLE_CHECKS
    data[ii] = numeric_cast<T>(Vec::get_entry(vec, ii));
#else
    data[ii] = Vec::get_entry(vec, ii);
#endif
  return data;
}


template <class RangeType, class SourceType>
typename std::enable_if<is_vector<SourceType>::value && is_vector<RangeType>::value, RangeType>::type
convert_to(const SourceType& source)
{
  const size_t size = source.size();
  auto ret = create<RangeType>(size);
  for (size_t ii = 0; ii < size; ++ii)
    ret[ii] =
#ifndef DXT_DISABLE_CHECKS
        numeric_cast<typename VectorAbstraction<RangeType>::S>(
#endif
            source[ii]
#ifndef DXT_DISABLE_CHECKS
#endif
        );
  return ret;
} // ... convert_to(...)


template <class V, class CharType, class CharTraits>
std::enable_if_t<Dune::XT::Common::is_vector<V>::value && !Dune::XT::Common::VectorAbstraction<V>::has_ostream,
                 std::basic_ostream<CharType, CharTraits>&>
operator<<(std::basic_ostream<CharType, CharTraits>& out, const V& vec)
{
  using Vector = Dune::XT::Common::VectorAbstraction<V>;
  if (vec.size() == 0)
    out << "[]";
  else if (vec.size() == 1)
    out << Vector::get_entry(vec, 0);
  else {
    out << "[" << Vector::get_entry(vec, 0);
    for (decltype(vec.size()) ii = 1; ii < vec.size(); ++ii)
      out << " " << Vector::get_entry(vec, ii);
    out << "]";
  }
  return out;
} // ... operator<<(...)


} // namespace Common
} // namespace XT
} // namespace Dune


namespace std {


/// clang 3.6 does not consider the overload in the ns for some reason during resultion of a call in gtest
template <class V, class Alloc, class CharType, class CharTraits>
std::basic_ostream<CharType, CharTraits>& operator<<(std::basic_ostream<CharType, CharTraits>& out,
                                                     const std::vector<V, Alloc>& vec)
{
  Dune::XT::Common::operator<<(out, vec);
  return out;
} // ... operator<<(...)


/// clang 3.6 does not consider the overload in the ns for some reason during resultion of a call in gtest
template <class CharType, class CharTraits, typename _Tp, std::size_t _Nm>
std::basic_ostream<CharType, CharTraits>& operator<<(std::basic_ostream<CharType, CharTraits>& out,
                                                     const std::array<_Tp, _Nm>& vec)
{
  Dune::operator<<(out, vec);
  return out;
} // ... operator<<(...)


} // namespace std

#endif // DUNE_XT_COMMON_VECTOR_HH
