// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2017 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_XT_GRID_TYPE_TRAITS_HH
#define DUNE_XT_GRID_TYPE_TRAITS_HH

#include <type_traits>

#include <dune/grid/common/entity.hh>
#include <dune/grid/common/intersection.hh>
#include <dune/grid/common/gridview.hh>


#if HAVE_DUNE_FEM
#include <dune/fem/gridpart/common/gridpart.hh>
#endif

#include <dune/xt/common/type_traits.hh>

#include <dune/xt/grid/grids.hh>

namespace Dune {
namespace XT {
namespace Grid {


template <class T>
struct is_intersection : public std::false_type
{
};

template <class G, class I>
struct is_intersection<Dune::Intersection<G, I>> : public std::true_type
{
  typedef std::remove_const_t<G> GridType;
};


template <class T>
struct is_grid : public std::false_type
{
};

template <>
struct is_grid<Dune::OneDGrid> : public std::true_type
{
};

template <int dim, class Coordinates>
struct is_grid<Dune::YaspGrid<dim, Coordinates>> : public std::true_type
{
};

#if HAVE_ALBERTA

template <int dim, int dimworld>
struct is_grid<Dune::AlbertaGrid<dim, dimworld>> : public std::true_type
{
};

#endif // HAVE_ALBERTA
#if HAVE_DUNE_ALUGRID

template <int dim, int dimworld, ALUGridElementType elType, ALUGridRefinementType refineType, class Comm>
struct is_grid<Dune::ALUGrid<dim, dimworld, elType, refineType, Comm>> : public std::true_type
{
};

#endif // HAVE_DUNE_ALUGRID
#if HAVE_DUNE_UGGRID || HAVE_UG

template <int dim>
struct is_grid<Dune::UGGrid<dim>> : public std::true_type
{
};

#endif // HAVE_DUNE_UGGRID || HAVE_UG

#if HAVE_DUNE_SPGRID
template <class ct, int dim, template <int> class Ref, class Comm>
struct is_grid<Dune::SPGrid<ct, dim, Ref, Comm>> : public std::true_type
{
};

#endif // HAVE_DUNE_SPGRID

template <class T>
struct is_grid_view : public std::false_type
{
};

template <class T>
struct is_grid_view<Dune::GridView<T>> : public std::true_type
{
};


namespace internal {


// the following did not work, so we need to use the messy approach below

// template <class T>
// struct is_grid_part : public std::false_type {};
// template <class T>
// struct is_grid_part<Dune::Fem::GridPartInterface<T>> : public std::true_type {};


template <class T>
struct is_grid_part_helper
{
  DXTC_has_typedef_initialize_once(Traits);
  static const bool is_candidate = DXTC_has_typedef(Traits)<T>::value;
};

template <class T>
struct is_subdomain_part_helper
{
  DXTC_has_typedef_initialize_once(GlobalGridPartType);
  static const bool value = DXTC_has_typedef(GlobalGridPartType)<T>::value;
};

} // namespace internal


template <class T, bool is_candidate = internal::is_grid_part_helper<T>::is_candidate>
struct is_grid_part : public std::false_type
{
};

#if HAVE_DUNE_FEM

template <class T>
struct is_grid_part<T, true> : public std::is_base_of<Dune::Fem::GridPartInterface<typename T::Traits>, T>
{
};


#endif // HAVE_DUNE_FEM


template <class T>
struct is_layer : public std::integral_constant<bool,
                                                is_grid_view<T>::value || is_grid_part<T>::value
                                                    || internal::is_subdomain_part_helper<T>::value>
{
};

template <class T>
struct is_alugrid : public std::false_type
{
};

template <class T>
struct is_conforming_alugrid : public std::false_type
{
};

#if HAVE_DUNE_ALUGRID

template <int dim, int dimworld, ALUGridElementType elType, ALUGridRefinementType refineType, class Comm>
struct is_alugrid<ALUGrid<dim, dimworld, elType, refineType, Comm>> : public std::true_type
{
};

template <int dim, int dimworld, ALUGridElementType elType, class Comm>
struct is_conforming_alugrid<ALUGrid<dim, dimworld, elType, Dune::conforming, Comm>> : public std::true_type
{
};

#endif // HAVE_DUNE_ALUGRID


template <class T,
          bool view = is_grid_view<T>::value,
          bool part = is_grid_part<T>::value,
          bool intersection = is_intersection<T>::value>
struct extract_grid : public std::false_type
{
};

template <class T>
struct extract_grid<T, true, false, false>
{
  typedef typename T::Grid type;
};

template <class T>
struct extract_grid<T, false, true, false>
{
  typedef typename T::GridType type;
};

template <class T>
struct extract_grid<T, false, false, true>
{
  typedef typename is_intersection<T>::GridType type;
};


template <class T>
using extract_grid_t = typename extract_grid<T>::type;


template <class T, bool view = is_grid_view<T>::value, bool part = is_grid_part<T>::value>
struct extract_intersection : public std::false_type
{
};

template <class T>
struct extract_intersection<T, true, false>
{
  typedef typename T::Intersection type;
};

template <class T>
struct extract_intersection<T, false, true>
{
  typedef typename T::IntersectionType type;
};


template <class T>
using extract_intersection_t = typename extract_intersection<T>::type;


template <class T, bool view = is_grid_view<T>::value, bool part = is_grid_part<T>::value>
struct extract_entity : public std::false_type
{
};

template <class T>
struct extract_entity<T, true, false>
{
  typedef typename T::template Codim<0>::Entity type;
};

template <class T>
struct extract_entity<T, false, true>
{
  typedef typename T::template Codim<0>::EntityType type;
};


template <class T>
using extract_entity_t = typename extract_entity<T>::type;


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_TYPE_TRAITS_HH
