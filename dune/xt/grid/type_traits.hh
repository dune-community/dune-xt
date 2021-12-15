// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2020)
//   René Fritze     (2016 - 2020)
//   Tobias Leibner  (2016 - 2020)

#ifndef DUNE_XT_GRID_TYPE_TRAITS_HH
#define DUNE_XT_GRID_TYPE_TRAITS_HH

#include <type_traits>

#include <dune/grid/common/entity.hh>
#include <dune/grid/common/intersection.hh>
#include <dune/grid/common/gridview.hh>

#include <dune/xt/common/type_traits.hh>

#if HAVE_ALBERTA
#  include <dune/xt/common/disable_warnings.hh>
#  include <dune/grid/albertagrid.hh>
#  include <dune/xt/common/reenable_warnings.hh>
#endif

#if HAVE_DUNE_ALUGRID
#  include <dune/alugrid/grid.hh>
#endif

#if HAVE_DUNE_SPGRID
#  include <dune/grid/spgrid.hh>
#  include <dune/grid/spgrid/dgfparser.hh>
#endif

#if HAVE_DUNE_UGGRID || HAVE_UG
#  include <dune/grid/uggrid.hh>
#endif

#include <dune/grid/onedgrid.hh>
#include <dune/grid/yaspgrid.hh>

namespace Dune {
namespace GridGlue {


// forward
template <typename P0, typename P1, int I, int O>
class Intersection;


} // namespace GridGlue
namespace XT::Grid {


namespace internal {


template <class T>
struct has_traits_helper
{
  DXTC_has_typedef_initialize_once(Traits);
  static constexpr bool is_candidate = DXTC_has_typedef(Traits)<T>::value;
};

// forward
template <class CouplingIntersectionType, class MacroIntersectionType>
class CouplingIntersectionWithCorrectNormal;

} // namespace internal


template <class T>
struct is_intersection : public std::false_type
{
  using GridType = std::false_type;
  using InsideElementType = std::false_type;
  using OutsideElementType = std::false_type;
};

template <class G, class I>
struct is_intersection<Dune::Intersection<G, I>> : public std::true_type
{
  using GridType = std::remove_const_t<G>;
  using InsideElementType = typename Dune::Intersection<G, I>::Entity;
  using OutsideElementType = typename Dune::Intersection<G, I>::Entity;
};

template <typename P0, typename P1, int I, int O>
struct is_intersection<Dune::GridGlue::Intersection<P0, P1, I, O>> : public std::true_type
{
  using GridType = typename Dune::GridGlue::Intersection<P0, P1, I, O>::InsideGridView::Grid;
  using InsideElementType = typename Dune::GridGlue::Intersection<P0, P1, I, O>::InsideEntity;
  using OutsideElementType = typename Dune::GridGlue::Intersection<P0, P1, I, O>::OutsideEntity;
};

template <class G, class I, typename P0, typename P1, int It, int O>
struct is_intersection<
    Dune::XT::Grid::internal::CouplingIntersectionWithCorrectNormal<Dune::GridGlue::Intersection<P0, P1, It, O>,
                                                                    Dune::Intersection<G, I>>> : public std::true_type
{
  using GridType = typename Dune::GridGlue::Intersection<P0, P1, It, O>::InsideGridView::Grid;
  using InsideElementType = typename Dune::GridGlue::Intersection<P0, P1, It, O>::InsideEntity;
  using OutsideElementType = typename Dune::GridGlue::Intersection<P0, P1, It, O>::OutsideEntity;
};


template <class T>
using extract_inside_element_t = typename is_intersection<T>::InsideElementType;

template <class T>
using extract_outside_element_t = typename is_intersection<T>::OutsideElementType;


template <class T>
struct is_grid : public std::false_type
{};

template <class T>
struct is_grid<const T> : public is_grid<T>
{};

template <>
struct is_grid<Dune::OneDGrid> : public std::true_type
{};

template <int dim, class Coordinates>
struct is_grid<Dune::YaspGrid<dim, Coordinates>> : public std::true_type
{};

#if HAVE_ALBERTA

template <int dim, int dimworld>
struct is_grid<Dune::AlbertaGrid<dim, dimworld>> : public std::true_type
{};

#endif // HAVE_ALBERTA
#if HAVE_DUNE_ALUGRID

template <int dim, int dimworld, ALUGridElementType elType, ALUGridRefinementType refineType, class Comm>
struct is_grid<Dune::ALUGrid<dim, dimworld, elType, refineType, Comm>> : public std::true_type
{};

template <int dim, int dimworld, ALU3dGridElementType elType, class Comm>
struct is_grid<Dune::ALU3dGrid<dim, dimworld, elType, Comm>> : public std::true_type
{};

#endif // HAVE_DUNE_ALUGRID
#if HAVE_DUNE_UGGRID || HAVE_UG

template <int dim>
struct is_grid<Dune::UGGrid<dim>> : public std::true_type
{};

#endif // HAVE_DUNE_UGGRID || HAVE_UG

#if HAVE_DUNE_SPGRID
template <class ct, int dim, template <int> class Ref, class Comm>
struct is_grid<Dune::SPGrid<ct, dim, Ref, Comm>> : public std::true_type
{};

#endif // HAVE_DUNE_SPGRID


template <class T, int codim = 0>
struct is_entity : public std::false_type
{};

template <int cd, int dim, class GridImp, template <int, int, class> class EntityImp>
struct is_entity<Dune::Entity<cd, dim, GridImp, EntityImp>, cd> : public std::true_type
{};


template <class T, bool candidate = internal::has_traits_helper<std::remove_const_t<T>>::is_candidate>
struct is_view : public std::false_type
{};

template <class T>
struct is_view<const T, true> : public is_view<T, true>
{};

template <class T>
struct is_view<T, true> : public std::is_base_of<Dune::GridView<typename T::Traits>, std::remove_const_t<T>>
{};

template <class T, bool is_candidate = internal::has_traits_helper<T>::is_candidate>
struct is_part : public std::false_type
{};


template <class T>
struct is_layer : public std::integral_constant<bool, is_view<T>::value || is_part<T>::value>
{};


template <class T>
struct is_yaspgrid : public std::false_type
{};

template <int dim, class Coordinates>
struct is_yaspgrid<YaspGrid<dim, Coordinates>> : public std::true_type
{};

template <class T>
struct is_uggrid : public std::false_type
{};

#if HAVE_DUNE_UGGRID || HAVE_UG

template <int dim>
struct is_uggrid<UGGrid<dim>> : public std::true_type
{};

#endif


template <class T>
struct is_alugrid : public std::false_type
{};

template <class T>
struct is_conforming_alugrid : public std::false_type
{};

template <class T>
struct is_simplex_alugrid : public std::false_type
{};

template <class T>
struct is_cube_alugrid : public std::false_type
{};


#if HAVE_DUNE_ALUGRID

template <int dim, int dimworld, ALUGridElementType elType, ALUGridRefinementType refineType, class Comm>
struct is_alugrid<ALUGrid<dim, dimworld, elType, refineType, Comm>> : public std::true_type
{};

template <int dim, int dimworld, ALUGridElementType elType, class Comm>
struct is_conforming_alugrid<ALUGrid<dim, dimworld, elType, Dune::conforming, Comm>> : public std::true_type
{};

template <int dim, int dimworld, ALUGridRefinementType refineType, class Comm>
struct is_simplex_alugrid<ALUGrid<dim, dimworld, ALUGridElementType::simplex, refineType, Comm>> : public std::true_type
{};

template <int dim, int dimworld, ALUGridRefinementType refineType, class Comm>
struct is_cube_alugrid<ALUGrid<dim, dimworld, ALUGridElementType::cube, refineType, Comm>> : public std::true_type
{};

#endif // HAVE_DUNE_ALUGRID


template <class T,
          bool view = is_view<T>::value,
          bool part = is_part<T>::value,
          bool intersection = is_intersection<T>::value,
          bool entity = is_entity<T>::value>
struct extract_grid : public AlwaysFalse<T>
{};

template <class T>
struct extract_grid<T, true, false, false, false>
{
  using type = std::decay_t<typename T::Grid>;
};

template <class T>
struct extract_grid<T, false, true, false, false>
{
  using type = std::decay_t<typename T::GridType>;
};

template <class T>
struct extract_grid<T, false, false, true, false>
{
  using type = std::decay_t<typename is_intersection<T>::GridType>;
};

template <int cd, int dim, class GridImp, template <int, int, class> class EntityImp>
struct extract_grid<Dune::Entity<cd, dim, GridImp, EntityImp>, false, false, false, true>
{
  using type = std::decay_t<GridImp>;
};

template <class T>
using extract_grid_t = typename extract_grid<T>::type;


template <class T, bool view = is_view<T>::value, bool part = is_part<T>::value>
struct extract_collective_communication : public AlwaysFalse<T>
{};

template <class T>
struct extract_collective_communication<T, true, false>
{
  using type = typename T::CollectiveCommunication;
};

template <class T>
struct extract_collective_communication<T, false, true>
{
  using type = typename T::CollectiveCommunicationType;
};

template <class T>
using extract_collective_communication_t = typename extract_collective_communication<T>::type;


template <class T, bool view = is_view<T>::value, bool part = is_part<T>::value>
struct extract_index_set : public AlwaysFalse<T>
{};

template <class T>
struct extract_index_set<T, true, false>
{
  using type = typename T::IndexSet;
};

template <class T>
struct extract_index_set<T, false, true>
{
  using type = typename T::IndexSetType;
};

template <class T>
using extract_index_set_t = typename extract_index_set<T>::type;


template <class T, bool view = is_view<T>::value, bool part = is_part<T>::value>
struct extract_intersection : public AlwaysFalse<T>
{};

template <class T>
struct extract_intersection<T, true, false>
{
  using type = typename T::Intersection;
};

template <class T>
struct extract_intersection<T, false, true>
{
  using type = typename T::IntersectionType;
};

template <class T>
using extract_intersection_t = typename extract_intersection<T>::type;


template <class T, bool view = is_view<T>::value, bool part = is_part<T>::value>
struct extract_intersection_iterator : public AlwaysFalse<T>
{};

template <class T>
struct extract_intersection_iterator<T, true, false>
{
  using type = typename T::IntersectionIterator;
};

template <class T>
struct extract_intersection_iterator<T, false, true>
{
  using type = typename T::IntersectionIteratorType;
};

template <class T>
using extract_intersection_iterator_t = typename extract_intersection_iterator<T>::type;


template <class T,
          size_t codim = 0,
          bool view = is_view<T>::value,
          bool part = is_part<T>::value,
          bool grid = is_grid<T>::value>
struct extract_entity : public AlwaysFalse<T>
{};

// template <class T, size_t codim, bool view, bool part>
// struct extract_entity<const T&, codim, view, part> : public extract_entity<T, codim, view, part>
//{};

template <class T, size_t codim>
struct extract_entity<T, codim, true, false, false>
{
  using type = typename T::template Codim<codim>::Entity;
};

template <class T, size_t codim>
struct extract_entity<T, codim, false, true, false>
{
  using type = typename T::template Codim<codim>::EntityType;
};

template <class T, size_t codim>
struct extract_entity<T, codim, false, false, true>
{
  using type = typename T::template Codim<codim>::Entity;
};

template <class T, size_t codim = 0>
using extract_entity_t = typename extract_entity<T, codim>::type;


template <class T, size_t codim = 0, bool view = is_view<T>::value, bool part = is_part<T>::value>
struct extract_local_geometry : public AlwaysFalse<T>
{};

template <class T, size_t codim>
struct extract_local_geometry<T, codim, true, false>
{
  using type = typename T::template Codim<codim>::LocalGeometry;
};

template <class T, size_t codim>
struct extract_local_geometry<T, codim, false, true>
{
  using type = typename T::template Codim<codim>::LocalGeometryType;
};

template <class T, size_t codim = 0>
using extract_local_geometry_t = typename extract_local_geometry<T, codim>::type;


template <class T, size_t codim = 0, bool view = is_view<T>::value, bool part = is_part<T>::value>
struct extract_geometry : public AlwaysFalse<T>
{};

template <class T, size_t codim>
struct extract_geometry<T, codim, true, false>
{
  using type = typename T::template Codim<codim>::Geometry;
};

template <class T, size_t codim>
struct extract_geometry<T, codim, false, true>
{
  using type = typename T::template Codim<codim>::GeometryType;
};

template <class T, size_t codim = 0>
using extract_geometry_t = typename extract_geometry<T, codim>::type;


template <class T,
          int c = 0,
          PartitionIteratorType pit = All_Partition,
          bool view = is_view<T>::value,
          bool part = is_part<T>::value>
struct extract_iterator : public AlwaysFalse<T>
{};

template <class T, int c, PartitionIteratorType pit>
struct extract_iterator<T, c, pit, true, false>
{
  using type = typename T::template Codim<c>::template Partition<pit>::Iterator;
};

template <class T, int c, PartitionIteratorType pit>
struct extract_iterator<T, c, pit, false, true>
{
  using type = typename T::template Codim<c>::template Partition<pit>::IteratorType;
};

template <class T, int c = 0, PartitionIteratorType pit = All_Partition>
using extract_iterator_t = typename extract_iterator<T, c, pit>::type;

//! struct to be used as comparison function e.g. in a std::map<Entity, EntityLess>
template <class GV>
struct EntityLess
{
  using IndexSet = typename GV::IndexSet;
  using E = typename GV::Grid::template Codim<0>::Entity;

  EntityLess(const IndexSet& index_set)
    : index_set_(index_set)
  {}

  bool operator()(const E& a, const E& b) const
  {
    return index_set_.index(a) < index_set_.index(b);
  }

  const IndexSet& index_set_;
};


} // namespace XT::Grid
} // namespace Dune

#endif // DUNE_XT_GRID_TYPE_TRAITS_HH
