// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2018)
//   Tobias Leibner  (2017 - 2018)

#ifndef DUNE_XT_LA_TYPE_TRAITS_HH
#define DUNE_XT_LA_TYPE_TRAITS_HH

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


static constexpr Backends default_backend = Backends::istl_sparse;
static constexpr Backends default_sparse_backend = Backends::istl_sparse;
static constexpr Backends default_dense_backend = Backends::common_dense;


// forwards
template <class Traits>
class ProvidesBackend;

template <class Traits>
class ProvidesDataAccess;

template <class Traits, class ScalarImp>
class ContainerInterface;

template <class Traits, class ScalarImp>
class MatrixInterface;

template <class Traits, class ScalarImp>
class VectorInterface;

template <class ScalarType, Backends backend>
struct Container;


namespace internal {


template <class C>
struct provides_backend_helper
{
  DXTC_has_typedef_initialize_once(Traits);

  static const bool is_candidate = DXTC_has_typedef(Traits)<C>::value;
}; // provides_backend_helper


template <class C>
struct provides_data_access_helper
{
  DXTC_has_typedef_initialize_once(Traits);
  DXTC_has_typedef_initialize_once(DataType);

  static const bool is_candidate = DXTC_has_typedef(Traits)<C>::value && DXTC_has_typedef(DataType)<C>::value;
}; // provides_data_access_helper


template <class C>
struct is_container_helper
{
  DXTC_has_typedef_initialize_once(Traits);
  DXTC_has_typedef_initialize_once(ScalarType);

  static const bool is_candidate = DXTC_has_typedef(Traits)<C>::value && DXTC_has_typedef(ScalarType)<C>::value;
}; // class is_container_helper


template <class M>
struct is_matrix_helper
{
  DXTC_has_typedef_initialize_once(Traits);
  DXTC_has_typedef_initialize_once(ScalarType);

  static const bool is_candidate = DXTC_has_typedef(Traits)<M>::value && DXTC_has_typedef(ScalarType)<M>::value;
}; // class is_matrix_helper


template <class V>
struct is_vector_helper
{
  DXTC_has_typedef_initialize_once(Traits);
  DXTC_has_typedef_initialize_once(ScalarType);

  static const bool is_candidate = DXTC_has_typedef(Traits)<V>::value && DXTC_has_typedef(ScalarType)<V>::value;
}; // class is_vector_helper


} // namespace internal


template <class C, bool candidate = internal::is_container_helper<C>::is_candidate>
struct is_container : public std::is_base_of<ContainerInterface<typename C::Traits, typename C::ScalarType>, C>
{};

template <class C>
struct is_container<C, false> : public std::false_type
{};


template <class M, bool candidate = internal::is_matrix_helper<M>::is_candidate>
struct is_matrix : public std::is_base_of<MatrixInterface<typename M::Traits, typename M::ScalarType>, M>
{};

template <class M>
struct is_matrix<M, false> : public std::false_type
{};


template <class V, bool candidate = internal::is_vector_helper<V>::is_candidate>
struct is_vector : public std::is_base_of<VectorInterface<typename V::Traits, typename V::ScalarType>, V>
{};

template <class V>
struct is_vector<V, false> : public std::false_type
{};


template <class C, bool candidate = internal::provides_backend_helper<C>::is_candidate>
struct provides_backend : public std::is_base_of<ProvidesBackend<typename C::Traits>, C>
{};

template <class C>
struct provides_backend<C, false> : public std::false_type
{};


template <class C, bool candidate = internal::provides_data_access_helper<C>::is_candidate>
struct provides_data_access : public std::is_base_of<ProvidesDataAccess<typename C::Traits>, C>
{};

template <class C>
struct provides_data_access<C, false> : public std::false_type
{};


template <class V, bool = is_vector<V>::value>
struct extract_matrix;

template <class V>
struct extract_matrix<V, false>
{
  static_assert(AlwaysFalse<V>::value, "V is not a vector!");
};

template <class V>
struct extract_matrix<V, true>
{
  using type = typename Container<typename V::ScalarType, V::sparse_matrix_type>::MatrixType;
};

template <class V>
using matrix_t = typename extract_matrix<V>::type;


template <class M, bool = is_matrix<M>::value>
struct extract_vector;

template <class M>
struct extract_vector<M, false>
{
  static_assert(AlwaysFalse<M>::value, "M is not a matrix!");
};

template <class M>
struct extract_vector<M, true>
{
  using type = typename Container<typename M::ScalarType, M::vector_type>::VectorType;
};

template <class M>
using vector_t = typename extract_vector<M>::type;


} // namespace LA
} // namespace XT
} // namespace Dune

#endif // #ifndef DUNE_XT_LA_TYPE_TRAITS_HH
