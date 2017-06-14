// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2017 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2017)
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


// forwards
template <class GlobalGridPartImp>
class SubdomainGridPart;

template <class GlobalGridPartImp>
class SubdomainBoundaryGridPart;

template <class GlobalGridPartImp>
class SubdomainCouplingGridPart;


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


template <class T, int codim = 0>
struct is_entity : public std::false_type
{
};

template <int cd, int dim, class GridImp, template <int, int, class> class EntityImp>
struct is_entity<Dune::Entity<cd, dim, GridImp, EntityImp>, cd> : public std::true_type
{
};


template <class T>
struct is_view : public std::false_type
{
};

template <class T>
struct is_view<Dune::GridView<T>> : public std::true_type
{
};

template <class T>
struct is_view<const T> : public is_view<T>
{
};


#if HAVE_DUNE_FEM

template <class T>
struct is_view<Dune::Fem::GridPart2GridView<T>> : public std::true_type
{
};

#endif // HAVE_DUNE_FEM

template <class T>
struct DUNE_DEPRECATED_MSG("Use is_view instead (03.04.2017)!") is_grid_view : public is_view<T>
{
};


namespace internal {


template <class T>
struct is_dd_subdomain_helper
{
  DXTC_has_typedef_initialize_once(GlobalGridPartType);
  static const bool value = DXTC_has_typedef(GlobalGridPartType)<T>::value;
};


// the following did not work, so we need to use the messy approach below

// template <class T>
// struct is_part : public std::false_type {};
// template <class T>
// struct is_part<Dune::Fem::GridPartInterface<T>> : public std::true_type {};

template <class T>
struct is_part_helper
{
  DXTC_has_typedef_initialize_once(Traits);
  static const bool is_candidate = DXTC_has_typedef(Traits)<T>::value;
};


} // namespace internal


template <class T, bool is_candidate = internal::is_dd_subdomain_helper<T>::value>
struct is_dd_subdomain : public std::false_type
{
};

template <class T>
struct is_dd_subdomain<T, true> : public std::is_base_of<SubdomainGridPart<typename T::GlobalGridPartType>, T>
{
};


template <class T, bool is_candidate = internal::is_dd_subdomain_helper<T>::value>
struct is_dd_subdomain_boundary : public std::false_type
{
};

template <class T>
struct is_dd_subdomain_boundary<T, true>
    : public std::is_base_of<SubdomainBoundaryGridPart<typename T::GlobalGridPartType>, T>
{
};


template <class T, bool is_candidate = internal::is_dd_subdomain_helper<T>::value>
struct is_dd_subdomain_coupling : public std::false_type
{
};

template <class T>
struct is_dd_subdomain_coupling<T, true>
    : public std::is_base_of<SubdomainCouplingGridPart<typename T::GlobalGridPartType>, T>
{
};


template <class T, bool is_candidate = internal::is_part_helper<T>::is_candidate>
struct is_part : public std::false_type
{
};

#if HAVE_DUNE_FEM

template <class T>
struct is_part<T, true> : public std::is_base_of<Dune::Fem::GridPartInterface<typename T::Traits>, T>
{
};

#endif // HAVE_DUNE_FEM

template <class T>
struct DUNE_DEPRECATED_MSG("Use is_part instead (03.04.2017)!") is_grid_part : public is_part<T>
{
};


template <class T>
struct is_layer : public std::integral_constant<bool,
                                                is_view<T>::value || is_part<T>::value || is_dd_subdomain<T>::value
                                                    || is_dd_subdomain_boundary<T>::value>
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
          bool view = is_view<T>::value,
          bool part = is_part<T>::value || is_dd_subdomain<T>::value || is_dd_subdomain_boundary<T>::value,
          bool intersection = is_intersection<T>::value>
struct extract_grid
{
  static_assert(Dune::AlwaysFalse<T>::value, "type T is not something we can extract a gridtype from");
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


template <class T,
          bool view = is_view<T>::value,
          bool part = is_part<T>::value || is_dd_subdomain<T>::value || is_dd_subdomain_boundary<T>::value>
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


template <class T,
          size_t codim = 0,
          bool view = is_view<T>::value,
          bool part = is_part<T>::value || is_dd_subdomain<T>::value || is_dd_subdomain_boundary<T>::value>
struct extract_entity : public std::false_type
{
};

template <class T, size_t codim>
struct extract_entity<T, codim, true, false>
{
  typedef typename T::template Codim<codim>::Entity type;
};

template <class T, size_t codim>
struct extract_entity<T, codim, false, true>
{
  typedef typename T::template Codim<codim>::EntityType type;
};


template <class T, size_t codim = 0>
using extract_entity_t = typename extract_entity<T, codim>::type;


template <class T,
          int c = 0,
          bool view = is_view<T>::value,
          bool part = is_part<T>::value || is_dd_subdomain<T>::value || is_dd_subdomain_boundary<T>::value>
struct extract_iterator : public std::false_type
{
};

template <class T, int c>
struct extract_iterator<T, c, true, false>
{
  typedef typename T::template Codim<c>::Iterator type;
};

template <class T, int c>
struct extract_iterator<T, c, false, true>
{
  typedef typename T::template Codim<c>::IteratorType type;
};


template <class T, int c = 0>
using extract_iterator_t = typename extract_iterator<T, c>::type;


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_TYPE_TRAITS_HH
