// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2017 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Tobias Leibner  (2017)

#ifndef DUNE_XT_LA_TYPE_TRAITS_HH
#define DUNE_XT_LA_TYPE_TRAITS_HH

#include <dune/xt/common/type_traits.hh>


namespace Dune {
namespace XT {
namespace LA {


// forwards
template <class Traits>
class ProvidesBackend;

template <class Traits>
class ProvidesDataAccess;

template <class Traits, class ScalarImp>
class ContainerInterface;


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


} // namespace internal


template <class C, bool candidate = internal::is_container_helper<C>::is_candidate>
struct is_container : public std::is_base_of<ContainerInterface<typename C::Traits, typename C::ScalarType>, C>
{
};

template <class C>
struct is_container<C, false> : public std::false_type
{
};


template <class C, bool candidate = internal::provides_backend_helper<C>::is_candidate>
struct provides_backend : public std::is_base_of<ProvidesBackend<typename C::Traits>, C>
{
};

template <class C>
struct provides_backend<C, false> : public std::false_type
{
};


template <class C, bool candidate = internal::provides_data_access_helper<C>::is_candidate>
struct provides_data_access : public std::is_base_of<ProvidesDataAccess<typename C::Traits>, C>
{
};

template <class C>
struct provides_data_access<C, false> : public std::false_type
{
};


} // namespace LA
} // namespace XT
} // namespace Dune

#endif // #ifndef DUNE_XT_LA_TYPE_TRAITS_HH
