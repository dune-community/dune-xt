// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2018)
//   René Fritze     (2015 - 2019)
//   Tobias Leibner  (2015 - 2018, 2020)

#ifndef DUNE_XT_GRID_VIEW_NEW_COUPLING_HH
#define DUNE_XT_GRID_VIEW_NEW_COUPLING_HH

#include <bitset>
#include <iterator>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <dune/geometry/typeindex.hh>

#include <dune/grid/common/gridview.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/grid/entity.hh>
#include <dune/xt/grid/search.hh>
#include <dune/xt/grid/type_traits.hh>

template <class type>
struct CompareType
{
  bool operator()(const type& one, const type& other) const
  {
    return (std::addressof(one) < std::addressof(other));
  }
};

namespace Dune::XT::Grid {
namespace internal {


template <class BaseGridViewType, bool codim_iters_provided, int codim>
class IndexMapCreator
{
public:
  using ElementType = XT::Grid::extract_entity_t<BaseGridViewType, 0>;
  using EntityType = XT::Grid::extract_entity_t<BaseGridViewType, codim>;
  using ElementPtrVectorType = std::vector<std::shared_ptr<const ElementType>>;
  using IndexType = typename extract_index_set_t<BaseGridViewType>::IndexType;
  static constexpr size_t dimDomain = BaseGridViewType::dimension;
  static constexpr size_t num_geometries = GlobalGeometryTypeIndex::size(dimDomain);
  using DomainType = typename BaseGridViewType::Intersection::GlobalCoordinate;

  IndexMapCreator(const DomainType& lower_left,
                  const DomainType& upper_right,
                  const std::bitset<dimDomain>& periodic_directions,
                  const BaseGridViewType& base_grid_view,
                  std::array<IndexType, dimDomain + 1>& entity_counts,
                  std::array<IndexType, num_geometries>& type_counts,
                  std::array<std::unordered_set<IndexType>, num_geometries>& entities_to_skip,
                  std::array<std::vector<IndexType>, num_geometries>& new_indices,
                  std::array<std::unordered_map<IndexType, ElementPtrVectorType>, num_geometries>&
                      boundary_entity_to_periodic_neighbors_maps)
    : lower_left_(lower_left)
    , upper_right_(upper_right)
    , periodic_directions_(periodic_directions)
    , base_grid_view_(base_grid_view)
    , base_index_set_(base_grid_view_.indexSet())
    , entity_counts_(entity_counts)
    , type_counts_(type_counts)
    , entities_to_skip_(entities_to_skip)
    , new_indices_(new_indices)
    , current_new_index_({})
    , boundary_entity_to_periodic_neighbors_maps_(boundary_entity_to_periodic_neighbors_maps)
  {
    for (const auto& geometry_type : base_index_set_.types(codim)) {
      const auto type_index = GlobalGeometryTypeIndex::index(geometry_type);
      const auto num_type_entities = base_index_set_.size(geometry_type);
      if constexpr (codim == 0)
        type_counts_[GlobalGeometryTypeIndex::index(geometry_type)] = num_type_entities;
      new_indices_[type_index].resize(num_type_entities);
    }
  }

  void create_index_map()
  {
    if constexpr (codim_iters_provided) {
      // if iterators for all codimensions are available, we can simply use entities(grid_view, codim)
      for (const auto& entity : entities(base_grid_view_, Dune::Codim<codim>())) {
        const auto index_in_base = base_grid_view_.indexSet().index(entity);
        const auto type_index = GlobalGeometryTypeIndex::index(entity.type());
        this->loop_body(entity, type_index, index_in_base);
      }
    } else {
      // In this case, there is only an iterator for the codim 0 entitys. To walk over the entities of codimension
      // codim, we thus walk over the respective subentities of the elements. As we are visiting each codim entity
      // several times this way, we have to remember the entities which we already visited and skip these when
      // encountering them again.
      for (const auto& element : elements(base_grid_view_)) {
        for (IndexType local_index = 0; local_index < element.subEntities(codim); ++local_index) {
          const auto& entity = element.template subEntity<codim>(local_index);
          const auto index_in_base = base_grid_view_.indexSet().index(entity);
          const auto type_index = GlobalGeometryTypeIndex::index(entity.type());
          if (!visited_entities_[type_index].count(index_in_base)) {
            this->loop_body(entity, type_index, index_in_base);
            visited_entities_[type_index].insert(index_in_base);
          } // if (entity has not been visited before)
        } // walk subentities in a given codimension
      } // walk codim0 elements
    }
    this->after_loop();
  } // ... create_index_map(...)

private:
  void loop_body(const EntityType& entity, const std::size_t& type_index, const IndexType& index_in_base)
  {
    if constexpr (codim == 0) {
      if (entity.hasBoundaryIntersections()) {
        // walk over intersections and store global coordinates of periodic neighbor entities (to be able to search
        // these entities later on via XT::Grid::EntitySearch)
        for (const auto& intersection : Dune::intersections(base_grid_view_, entity)) {
          const int index_in_inside = intersection.indexInInside();
          if (intersection.boundary()) {
            bool is_periodic = false;
            auto periodic_neighbor_coords = intersection.geometry().center();
            size_t num_boundary_coords = 0;
            for (std::size_t ii = 0; ii < dimDomain; ++ii) {
              if (periodic_directions_[ii]) {
                if (XT::Common::FloatCmp::eq(periodic_neighbor_coords[ii], lower_left_[ii])) {
                  is_periodic = true;
                  periodic_neighbor_coords[ii] =
                      upper_right_[ii] - 1.0 / 100.0 * (entity.geometry().center()[ii] - lower_left_[ii]);
                  ++num_boundary_coords;
                } else if (XT::Common::FloatCmp::eq(periodic_neighbor_coords[ii], upper_right_[ii])) {
                  is_periodic = true;
                  periodic_neighbor_coords[ii] =
                      lower_left_[ii] + 1.0 / 100.0 * (upper_right_[ii] - entity.geometry().center()[ii]);
                  ++num_boundary_coords;
                }
              }
            }
            if (is_periodic) {
              assert(num_boundary_coords == 1);
              periodic_coords_.push_back(periodic_neighbor_coords);
              periodic_coords_index_.push_back(std::make_tuple(type_index, index_in_base, index_in_inside));
            }
          }
        }
        boundary_entity_to_periodic_neighbors_maps_[type_index][index_in_base] =
            ElementPtrVectorType(entity.subEntities(1), nullptr);
      } // if (entity.hasBoundaryIntersections)
    } else {
      // check if entity is on right-hand periodic boundary, in that case it will be identified with its periodically
      // equivalent entity on the left
      auto periodic_coords = entity.geometry().center();
      std::size_t num_upper_right_coords = 0;
      for (std::size_t ii = 0; ii < dimDomain; ++ii) {
        if (periodic_directions_[ii]) {
          if (XT::Common::FloatCmp::eq(periodic_coords[ii], upper_right_[ii])) {
            ++num_upper_right_coords;
            periodic_coords[ii] = lower_left_[ii];
          }
        }
      }

      if (num_upper_right_coords == 0) {
        // assign new index and increase index counter
        new_indices_[type_index][index_in_base] = current_new_index_[type_index];
        ++current_new_index_[type_index];
        // increase counters for geometry type and entities of this codimension
        ++type_counts_[type_index];
        ++entity_counts_[codim];
      } else {
        // entity is on right periodic boundary and will be ignored
        entities_to_skip_[type_index].insert(index_in_base);
        // store global coordinates of this entity to be able to ignore it later on
        periodic_coords_.push_back(periodic_coords);
        periodic_coords_index_.push_back({type_index, index_in_base, 0});
      }
    }
  } // loop_body

  void after_loop()
  {
    if constexpr (codim == 0)
      entity_counts_[codim] = base_index_set_.size(0);

    // now find periodic neighbor entities (for the intersections, if codim == 0) or periodically equivalent entities
    // (codim != 0)
    typename std::conditional<codim_iters_provided,
                              EntityInlevelSearch<BaseGridViewType, codim>,
                              FallbackEntityInlevelSearch<BaseGridViewType, codim>>::type
        entity_search_codim(base_grid_view_);
    auto periodic_entity_ptrs = entity_search_codim(periodic_coords_);
    for (size_t vector_index = 0; vector_index < periodic_entity_ptrs.size(); ++vector_index) {
      const auto& index = periodic_coords_index_[vector_index];
      auto& periodic_entity_ptr = periodic_entity_ptrs[vector_index];
      if (!periodic_entity_ptr)
        DUNE_THROW(Dune::InvalidStateException, "Could not find periodic neighbor entity");
      const auto& type_index = std::get<0>(index);
      const auto& entity_index = std::get<1>(index);
      if constexpr (codim == 0) {
        // store periodic neighbor for intersection
        const auto& local_intersection_index = std::get<2>(index);
        boundary_entity_to_periodic_neighbors_maps_[type_index][entity_index][local_intersection_index] =
            std::move(periodic_entity_ptr);
      } else {
        // assign index of periodic equivalent entity to entities that are replaced
        const auto periodic_entity_index = base_index_set_.index(*periodic_entity_ptr);
        const auto& periodic_entity_type_index = GlobalGeometryTypeIndex::index(periodic_entity_ptr->type());
        new_indices_[type_index][entity_index] = new_indices_[periodic_entity_type_index][periodic_entity_index];
      }
    }
  } // after_loop()

  const DomainType& lower_left_;
  const DomainType& upper_right_;
  const std::bitset<dimDomain>& periodic_directions_;
  const BaseGridViewType& base_grid_view_;
  const extract_index_set_t<BaseGridViewType>& base_index_set_;
  std::array<IndexType, dimDomain + 1>& entity_counts_;
  std::array<IndexType, num_geometries>& type_counts_;
  std::array<std::unordered_set<IndexType>, num_geometries>& entities_to_skip_;
  std::array<std::vector<IndexType>, num_geometries>& new_indices_;
  std::vector<DomainType> periodic_coords_;
  std::vector<std::tuple<size_t, IndexType, int>> periodic_coords_index_;
  std::array<IndexType, num_geometries> current_new_index_;
  std::array<std::unordered_map<IndexType, ElementPtrVectorType>, num_geometries>&
      boundary_entity_to_periodic_neighbors_maps_;
  std::array<std::unordered_set<IndexType>, GlobalGeometryTypeIndex::size(dimDomain)> visited_entities_;
}; // struct IndexMapCreator< ... >

/** \brief IndexSet for CouplingGridViewWrapper
 *
 * PeriodicIndexSet is derived from the IndexSet of the underlying grid view. Other than the IndexSet,
 * PeriodicIndexSet returns the same index for entities that are periodically equivalent, i.e. entities on periodic
 * boundaries that are regarded as the same entity in the periodic setting. Consequently, the PeriodicIndexSet is
 * usually smaller than the IndexSet and the size(...) methods return lower values than the corresponding methods of
 * the non-periodic IndexSet.
 *
 * \see CouplingGridView
 */
template <class BaseGridViewImp>
class PeriodicIndexSet
  : public Dune::IndexSet<extract_grid_t<BaseGridViewImp>,
                          PeriodicIndexSet<BaseGridViewImp>,
                          typename extract_index_set_t<BaseGridViewImp>::IndexType,
                          typename extract_index_set_t<BaseGridViewImp>::Types>
{
  using BaseGridViewType = BaseGridViewImp;
  using ThisType = PeriodicIndexSet;
  using BaseIndexSetType = extract_index_set_t<BaseGridViewType>;
  using BaseType = Dune::IndexSet<extract_grid_t<BaseGridViewType>,
                                  ThisType,
                                  typename BaseIndexSetType::IndexType,
                                  typename BaseIndexSetType::Types>;

public:
  using BaseType::dimension;
  using typename BaseType::IndexType;
  using typename BaseType::Types;
  static constexpr size_t num_geometries = GlobalGeometryTypeIndex::size(dimension);

  PeriodicIndexSet(const BaseIndexSetType& base_index_set,
                   const std::array<IndexType, dimension + 1>& entity_counts,
                   const std::array<IndexType, num_geometries>& type_counts,
                   const std::array<std::vector<IndexType>, num_geometries>& new_indices)
    : BaseType()
    , base_index_set_(base_index_set)
    , entity_counts_(entity_counts)
    , type_counts_(type_counts)
    , new_indices_(new_indices)
  {
    assert(entity_counts_.size() >= dimension + 1);
  }

  PeriodicIndexSet(const ThisType& other) = default;
  PeriodicIndexSet(ThisType&& other) = default;
  // assigment currently does not work due to the reference members
  ThisType& operator=(const ThisType& other) = delete;
  ThisType& operator=(ThisType&& other) = delete;


  template <int cd, class CodimCdEntityType>
  IndexType index(const CodimCdEntityType& entity) const
  {
    IndexType base_entity_index = base_index_set_.template index<cd>(entity);
    if constexpr (cd == 0)
      return base_entity_index;
    else {
      const auto type_index = GlobalGeometryTypeIndex::index(entity.type());
      return new_indices_[type_index][base_entity_index];
    }
  }

  template <class EntityType>
  IndexType index(const EntityType& entity) const
  {
    return index<EntityType::codimension>(entity);
  }

  template <int cd, class CodimCdEntityType>
  IndexType subIndex(const CodimCdEntityType& entity, int i, unsigned int codim) const
  {
    IndexType base_sub_index = base_index_set_.template subIndex<cd>(entity, i, codim);
    if (codim == 0)
      return base_sub_index;
    else {
      const auto& ref_element = reference_element(entity);
      const auto type_index = GlobalGeometryTypeIndex::index(ref_element.type(i, codim));
      return new_indices_[type_index][base_sub_index];
    }
  }

  template <class EntityType>
  IndexType subIndex(const EntityType& entity, int i, unsigned int codim) const
  {
    return subIndex<EntityType::codimension>(entity, i, codim);
  }

  Types types(int codim) const
  {
    return base_index_set_.types(codim);
  }

  IndexType size(Dune::GeometryType type) const
  {
    const auto type_index = GlobalGeometryTypeIndex::index(type);
    return type_counts_[type_index];
  }

  IndexType size(int codim) const
  {
    assert(codim <= dimension);
    return entity_counts_[codim];
  }

  template <class EntityType>
  bool contains(const EntityType& entity) const
  {
    return base_index_set_.contains(entity);
  }

private:
  const BaseIndexSetType& base_index_set_;
  const std::array<IndexType, dimension + 1>& entity_counts_;
  const std::array<IndexType, num_geometries>& type_counts_;
  const std::array<std::vector<IndexType>, num_geometries>& new_indices_;
}; // class PeriodicIndexSet<...>


/** \brief Intersection implementation for CouplingGridViewWrapper
 *
 * PeriodicIntersectionImp is derived from the Intersection of the underlying grid view. On the inside of the grid or
 * if periodic_ is false, the PeriodicIntersection will behave exactly like its BaseType. If periodic_ is true, the
 * PeriodicIntersection will return neighbor == true even if it actually is on the boundary. In this case, outside(),
 * geometryInOutside() and indexInOutside() are well-defined and give the information from the periodically adjacent
 * entity.
 *
 * \see CouplingGridView
 */
template <class BaseIntersectionImp>
class PeriodicIntersectionImp : public BaseIntersectionImp
{
  using BaseType = BaseIntersectionImp;

public:
  using typename BaseType::LocalGeometry;
  using ElementType = typename BaseType::Entity;

  //! \brief Constructor from base intersection
  template <class BaseGridViewType>
  PeriodicIntersectionImp(BaseType base_intersection,
                          const BaseGridViewType& base_grid_view,
                          const ElementType* periodic_neighbor)
    : BaseType(base_intersection)
    , periodic_neighbor_(periodic_neighbor)
  {
    if (periodic_neighbor_)
      intersection_in_outside_ = find_intersection_in_outside(base_grid_view);
  }

  //! \brief Default constructor
  PeriodicIntersectionImp()
    : BaseType()
  {}

  PeriodicIntersectionImp(const PeriodicIntersectionImp& other) = default;
  PeriodicIntersectionImp(PeriodicIntersectionImp&& other) = default;
  PeriodicIntersectionImp& operator=(const PeriodicIntersectionImp& other) = default;
  PeriodicIntersectionImp& operator=(PeriodicIntersectionImp&& other) = default;

  // methods that differ from BaseType
  bool neighbor() const
  {
    if (periodic_neighbor_)
      return true;
    else
      return BaseType::neighbor();
  } // bool neighbor() const

  ElementType outside() const
  {
    if (periodic_neighbor_)
      return *periodic_neighbor_;
    else
      return ElementType(BaseType::outside());
  } // ... outside() const

  LocalGeometry geometryInOutside() const
  {
    if (periodic_neighbor_)
      return intersection_in_outside_->geometryInInside();
    else
      return BaseType::geometryInOutside();
  } // ... geometryInOutside() const

  int indexInOutside() const
  {
    if (periodic_neighbor_)
      return intersection_in_outside_->indexInInside();
    else
      return BaseType::indexInOutside();
  } // int indexInOutside() const

private:
  // tries to find intersection in outside (works only if this intersection is on the periodic boundary)
  template <class BaseGridViewType>
  std::shared_ptr<const BaseType> find_intersection_in_outside(const BaseGridViewType& base_grid_view) const
  {
    const auto coords = this->geometry().center();
    // walk over outside intersections and find an intersection on the boundary that differs only in one coordinate
    for (const BaseType& other_intersection : Dune::intersections(base_grid_view, *periodic_neighbor_)) {
      if (other_intersection.boundary()) {
        const auto coord_diff = other_intersection.geometry().center() - coords;
        size_t coord_diff_count = 0;
        for (const auto& entry : coord_diff)
          if (Common::FloatCmp::ne(entry, 0.))
            ++coord_diff_count;
        if (coord_diff_count == 1)
          return std::make_shared<const BaseType>(other_intersection);
      }
    }
    DUNE_THROW(Dune::InvalidStateException, "Could not find outside intersection!");
    return nullptr;
  } // ... find_intersection_in_outside() const

protected:
  const ElementType* periodic_neighbor_;
  std::shared_ptr<const BaseType> intersection_in_outside_;
}; // ... class PeriodicIntersectionImp ...


/** \brief IntersectionIterator for CouplingGridViewWrapper
 *
 * PeriodicIntersectionIterator is derived from the IntersectionIterator of the underlying grid view and behaves
 * exactly like the underlying IntersectionIterator except that it returns a PeriodicIntersection in its operator* and
 * operator-> methods.
 *
 * \see CouplingGridView
 */
template <class BaseGridViewImp>
class PeriodicIntersectionIterator : public extract_intersection_iterator_t<BaseGridViewImp>
{
  using BaseGridViewType = BaseGridViewImp;
  using BaseType = extract_intersection_iterator_t<BaseGridViewImp>;
  using ThisType = PeriodicIntersectionIterator;
  using IntersectionImp = PeriodicIntersectionImp<extract_intersection_t<BaseGridViewImp>>;

public:
  using BaseIntersectionType = typename BaseType::Intersection;
  using Intersection = Dune::Intersection<extract_grid_t<BaseGridViewImp>, IntersectionImp>;
  using ElementType = extract_entity_t<BaseGridViewType, 0>;
  using ElementPtrVectorType = std::vector<std::shared_ptr<const ElementType>>;

  PeriodicIntersectionIterator(BaseType base_intersection_iterator,
                               const BaseGridViewType& base_grid_view,
                               const ElementType& element,
                               const ElementPtrVectorType* periodic_neighbors)
    : BaseType(base_intersection_iterator)
    , base_grid_view_(base_grid_view)
    , element_(element)
    , has_boundary_intersections_(element_.hasBoundaryIntersections())
    , periodic_neighbors_(periodic_neighbors)
    , current_intersection_(nullptr)
  {}

  PeriodicIntersectionIterator(const ThisType& other)
    : BaseType(other)
    , base_grid_view_(other.base_grid_view_)
    , element_(other.element_)
    , has_boundary_intersections_(other.has_boundary_intersections_)
    , periodic_neighbors_(other.periodic_neighbors_)
    , current_intersection_(nullptr)
  {}

  ThisType& operator=(const ThisType& other)
  {
    BaseType::operator=(other);
    base_grid_view_ = other.base_grid_view_;
    element_ = other.element_;
    has_boundary_intersections_ = other.has_boundary_intersections_;
    periodic_neighbors_ = other.periodic_neighbors_;
    current_intersection_ = nullptr;
    return *this;
  }

  PeriodicIntersectionIterator(ThisType&& other) = default;
  ThisType& operator=(ThisType&& other) = default;

  // methods that differ from BaseType
  Intersection operator*() const
  {
    return *create_current_intersection();
  }

  const Intersection* operator->() const
  {
    current_intersection_ = create_current_intersection();
    return current_intersection_.get();
  }

private:
  // create current intersection without checking if this is the end iterator
  auto create_current_intersection() const
  {
    const auto base_intersection = BaseType::operator*();
    const int index = base_intersection.indexInInside();
    assert(!has_boundary_intersections_ || periodic_neighbors_->size() > static_cast<size_t>(index));
    return std::make_unique<const Intersection>(
        IntersectionImp(base_intersection,
                        base_grid_view_,
                        has_boundary_intersections_ ? (*periodic_neighbors_)[index].get() : nullptr));
  } // ... create_current_intersection() const

  const BaseGridViewType& base_grid_view_;
  const ElementType& element_;
  const bool has_boundary_intersections_;
  const ElementPtrVectorType* periodic_neighbors_;
  mutable std::unique_ptr<const Intersection> current_intersection_;
}; // ... class PeriodicIntersectionIterator ...

///////// START HERE //////////

// forward for Traits
template <class BaseGridViewImp, class GridGlueType>
class CouplingGridViewWrapper;

//! Traits for CouplingGridViewWrapper
template <class BaseGridViewImp, class GridGlueType>
class CouplingGridViewWrapperTraits : public BaseGridViewImp::Traits
{
public:
  using BaseGridViewType = BaseGridViewImp;
  // use types from BaseGridViewType...
  using GridViewImp = CouplingGridViewWrapper<BaseGridViewType, GridGlueType>;
  using Grid = extract_grid_t<BaseGridViewType>;
  using IndexSet = PeriodicIndexSet<BaseGridViewType>;
  using CollectiveCommunication = extract_collective_communication_t<BaseGridViewType>;
  using BaseGridViewTraits = typename BaseGridViewType::Traits;
  using LocalElementType = typename GridGlueType::MicroEntityType;
//  static constexpr size_t num_geometries = GlobalGeometryTypeIndex::size(BaseGridViewImp::dimension);

  using GlueType = typename GridGlueType::GlueType;
  using CouplingIntersectionType = typename GlueType::Intersection;
  using Intersection = CouplingIntersectionType;
  using IntersectionIterator = typename std::set<CouplingIntersectionType, CompareType<CouplingIntersectionType>>::const_iterator;

  template <int cd>
  struct Codim : public BaseGridViewTraits::template Codim<cd>
  {
    // We need to define these in case BaseGridViewImp is a grid part.
    using Entity = extract_entity_t<BaseGridViewImp, cd>;
    using Geometry = extract_geometry_t<BaseGridViewImp, cd>;
    using LocalGeometry = extract_local_geometry_t<BaseGridViewImp, cd>;

    using Iterator = typename std::vector<LocalElementType>::const_iterator;

//    /* PeriodicIterator is the same as the Iterator of the BaseGridViewType, except that it visits only one entity of
//     * several periodically equivalent entities */
//    template <PartitionIteratorType pit>
//    class PeriodicIterator : public extract_iterator_t<BaseGridViewType, cd, pit>
//    {
//      using BaseType = extract_iterator_t<BaseGridViewType, cd, pit>;
//      using BaseIndexSetType = extract_index_set_t<BaseGridViewType>;

//    public:
//      using IndexType = typename IndexSet::IndexType;
//      using difference_type = std::ptrdiff_t;
//      using value_type = const extract_entity_t<BaseGridViewImp, cd>;
//      using pointer = value_type*;
//      using reference = value_type&;
//      using iterator_category = std::forward_iterator_tag;

//      PeriodicIterator(BaseType base_iterator,
//                       const std::array<std::unordered_set<IndexType>, num_geometries>* entities_to_skip,
//                       const BaseIndexSetType* base_index_set,
//                       const BaseType& base_it_end)
//        : BaseType(base_iterator)
//        , entities_to_skip_(entities_to_skip)
//        , base_index_set_(base_index_set)
//        , base_it_end_(std::make_shared<const BaseType>(base_it_end))
//      {}

//      PeriodicIterator& operator++()
//      {
//        BaseType::operator++();
//        while (cd > 0 && *this != *base_it_end_
//               && (*entities_to_skip_)[GlobalGeometryTypeIndex::index((*this)->type())].count(
//                   base_index_set_->index(this->operator*())))
//          BaseType::operator++();
//        return *this;
//      }

//      PeriodicIterator operator++(int)
//      {
//        return this->operator++();
//      }

//    private:
//      const std::array<std::unordered_set<IndexType>, num_geometries>* entities_to_skip_;
//      const BaseIndexSetType* base_index_set_;
//      std::shared_ptr<const BaseType> base_it_end_;
//    };

//    using Iterator = PeriodicIterator<All_Partition>;

    template <PartitionIteratorType pit>
    struct Partition : public BaseGridViewTraits::template Codim<cd>::template Partition<pit>
    {
      using Iterator = typename std::vector<LocalElementType>::const_iterator;
    }; // struct Partition
  }; // ... struct Codim ...

  enum
  {
    conforming = BaseGridViewTraits::conforming
  };
  enum
  {
    dimension = Grid::dimension
  };
  enum
  {
    dimensionworld = Grid::dimensionworld
  };

  using ctype = typename Grid::ctype;

  // ...except for the Intersection and IntersectionIterator
//  using Intersection = Dune::Intersection<Grid, PeriodicIntersectionImp<extract_intersection_t<BaseGridViewType>>>;
//  using IntersectionIterator = PeriodicIntersectionIterator<BaseGridViewType>;
}; // ... class CouplingGridViewWrapperTraits ...


/** \brief Actual Implementation of CouplingGridView
 *  \see CouplingGridView
 */
template <class BaseGridViewImp, class GridGlueType>
class CouplingGridViewWrapper : public BaseGridViewImp
{
  using BaseType = BaseGridViewImp;
  using ThisType = CouplingGridViewWrapper;
  using Traits = CouplingGridViewWrapperTraits<BaseType, GridGlueType>;

public:
  using GlueType = typename GridGlueType::GlueType;
  using CouplingIntersectionType = typename GlueType::Intersection;
  using IntersectionType = CouplingIntersectionType;

  using Grid = extract_grid_t<BaseType>;
  using IndexSet = PeriodicIndexSet<BaseType>;
  using ElementType = extract_entity_t<BaseType, 0>;
//  using IntersectionIterator = PeriodicIntersectionIterator<BaseType>;
//  using BaseIntersectionType = typename IntersectionIterator::BaseIntersectionType;
  using IndexType = typename Traits::IndexSet::IndexType;
//  using DomainType = typename BaseIntersectionType::GlobalCoordinate;
//  using Intersection = Dune::Intersection<Grid, PeriodicIntersectionImp<BaseIntersectionType>>;
//  using ElementPtrVectorType = std::vector<std::shared_ptr<const ElementType>>;
//  static constexpr size_t dimDomain = BaseType::dimension;

  using MacroGridViewType = typename GridGlueType::MacroGridViewType;
  using MacroGridType = typename GridGlueType::MacroGridType;
  using MacroElementType = typename MacroGridType::template Codim<0>::Entity;
  using MacroIntersectionType = typename MacroGridViewType::Intersection;

  using LocalGridProviderType = typename GridGlueType::LocalGridProviderType;
  using LocalGridType = typename GridGlueType::LocalGridType;
  using LocalGridViewType = typename GridGlueType::MicroGridViewType;
  using LocalElementType = typename GridGlueType::MicroEntityType;

  // TODO: add the macro intersection to use CorrectedCouplingIntersection
//    using CorrectedCouplingIntersectionType = CouplingIntersectionWithCorrectNormal<CouplingIntersectionType, MacroIntersectionType>;

  template <int cd>
  struct Codim : public Traits::template Codim<cd>
  {};

//private:
//  // compile time for loop to loop over the codimensions in constructor, see http://stackoverflow.com/a/11081785
//  template <int codim, int to>
//  struct static_for_loop_for_index_maps
//  {
//    template <class... Args>
//    void operator()(Args&&... args)
//    {
//      if constexpr (codim < to) {
//        IndexMapCreator<BaseType, codim_iters_provided, codim> index_map_creator(std::forward<Args>(args)...);
//        index_map_creator.create_index_map();
//        static_for_loop_for_index_maps<codim + 1, to>()(std::forward<Args>(args)...);
//      }
//    }
//  };

public:
  CouplingGridViewWrapper(const BaseType& base_grid_view,
                          const MacroElementType& ss,
                          const MacroElementType& nn,
                          GridGlueType& dd_grid)
    : BaseType(base_grid_view),
      inside_element_(ss),
      outside_element_(nn),
      dd_grid_(dd_grid),
      macro_grid_view_(dd_grid.macro_grid_view()),
      local_inside_grid_(dd_grid.local_grid(ss))

  {
    if (base_grid_view.comm().size() > 1)
      DUNE_THROW(Dune::NotImplemented, "CouplingGridView does not work for MPI-distributed grids!");
    this->update();
  } // constructor CouplingGridViewWrapper(...)

  CouplingGridViewWrapper(const ThisType& other) = default;
  CouplingGridViewWrapper(ThisType&& other) = default;
  ThisType& operator=(const ThisType& other) = default;
  ThisType& operator=(ThisType&& other) = default;

  // This is the only member function that actually changes the state of the CouplingGridView. The copy constructor and
  // copy assignment operators only do a shallow copy of the shared_ptrs, so we have to make sure that we reassign the
  // shared_ptrs in this function to avoid changing state of a copied-from CouplingGridView.
  void update()
  {
      auto& coupling = dd_grid_.coupling(inside_element_, -1, outside_element_, -1, true);
      for (auto coupling_intersection_it = coupling.template ibegin<0>();
           coupling_intersection_it != coupling.template iend<0>();
           ++coupling_intersection_it) {
          auto inside = coupling_intersection_it->inside();
          inside_elements_ids_.push_back(local_inside_grid_->leaf_view().indexSet().index(inside));
      }
      // some coupling intersection may have the same inside element, remove duplicates
      std::sort(inside_elements_ids_.begin(), inside_elements_ids_.end());
      auto last = std::unique(inside_elements_ids_.begin(), inside_elements_ids_.end());
      inside_elements_ids_.erase(last, inside_elements_ids_.end());

      for (auto id = inside_elements_ids_.begin(); id != inside_elements_ids_.end(); ++id) {
          for (auto el = local_inside_grid_->leaf_view().template begin<0>();
               el != local_inside_grid_->leaf_view().template end<0>(); ++el) {
              if (local_inside_grid_->leaf_view().indexSet().index(*el) == *id) {
                  // This is the inside element we are searching for.. add it to the vector
                  inside_elements_.push_back(*el);
                  std::set<CouplingIntersectionType, CompareType<CouplingIntersectionType>> coupling_intersection_set;
                  // now iteratate over all intersections to find all coupling intersections
                  for (auto coupling_intersection_it = coupling.template ibegin<0>();
                       coupling_intersection_it != coupling.template iend<0>();
                       ++coupling_intersection_it) {
                      auto inside = coupling_intersection_it->inside();
                      auto inside_id = local_inside_grid_->leaf_view().indexSet().index(inside);
                      if (inside_id == *id) {
//                            CorrectedCouplingIntersectionType coupling_intersection(*coupling_intersection_it, macro_intersection_);
                          coupling_intersection_set.insert(*coupling_intersection_it);
                      }
                  }
                  coupling_intersections_.push_back(coupling_intersection_set);
              }
          }
      }
      // introduce a local to global map
      for (auto id = inside_elements_ids_.begin(); id != inside_elements_ids_.end(); ++id) {
          for (auto el = local_inside_grid_->leaf_view().template begin<0>();
               el != local_inside_grid_->leaf_view().template end<0>(); ++el) {
              if (local_inside_grid_->leaf_view().indexSet().index(*el) == *id) {
                  // This is the inside element we are searching for..
                  local_to_inside_indices_.push_back({*id, local_to_inside_indices_.size()});
              }
          }
      }
  }

  size_t local_to_inside_index(const LocalElementType& local_element) const
  {
      auto id = local_inside_grid_->leaf_view().indexSet().index(local_element);
      for (auto index_pair = local_to_inside_indices_.begin() ; index_pair != local_to_inside_indices_.end(); ++index_pair) {
          if (id == index_pair->first) {
              return index_pair->second;
          }
      }
  }

  int size(int codim) const
  {
    DUNE_THROW(NotImplemented,"");
  }

  int size(const Dune::GeometryType& type) const
  {
    DUNE_THROW(NotImplemented,"");
  }

//  template<int cd>  // fake template
//  typename std::vector<LocalElementType>::iterator begin() { return inside_elements_.begin(); };

//  template<int cd>  // fake template
//  typename std::vector<LocalElementType>::iterator end() { return inside_elements_.end(); };

  typename std::set<CouplingIntersectionType, CompareType<CouplingIntersectionType>>::const_iterator ibegin(const LocalElementType& inside_element) const
  {
      return coupling_intersections_[local_to_inside_index(inside_element)].begin();
  };
  typename std::set<CouplingIntersectionType, CompareType<CouplingIntersectionType>>::const_iterator iend(const LocalElementType& inside_element) const
  {
      return coupling_intersections_[local_to_inside_index(inside_element)].end();
  };

  template <int cd>
  typename Codim<cd>::Iterator begin() const
  {
      return inside_elements_.begin();
  }

  template <int cd>
  typename Codim<cd>::Iterator end() const
  {
      return inside_elements_.end();
  }

  template <int cd, PartitionIteratorType pitype>
  typename Codim<cd>::template Partition<pitype>::Iterator begin() const
  {
    return inside_elements_.begin();
//    return typename Codim<cd>::template Partition<pitype>::Iterator(BaseType::template begin<cd, pitype>(),
//                                                                    entities_to_skip_.get(),
//                                                                    &BaseType::indexSet(),
//                                                                    BaseType::template end<cd, pitype>());
  }


  template <int cd, PartitionIteratorType pitype>
  typename Codim<cd>::template Partition<pitype>::Iterator end() const
  {
    return inside_elements_.end();
//    return typename Codim<cd>::template Partition<pitype>::Iterator(BaseType::template end<cd, pitype>(),
//                                                                    entities_to_skip_.get(),
//                                                                    &BaseType::indexSet(),
//                                                                    BaseType::template end<cd, pitype>());
  }

  const IndexSet& indexSet() const
  {
    DUNE_THROW(NotImplemented,"");
  }

//  IntersectionIterator ibegin(const typename Codim<0>::Entity& entity) const
//  {
//    const auto& type_index = GlobalGeometryTypeIndex::index(entity.type());
//    assert(!entity.hasBoundaryIntersections()
//           || (*boundary_entity_to_periodic_neighbors_maps_)[type_index].count(index_set_->index(entity)));
//    return IntersectionIterator(BaseType::ibegin(entity),
//                                *this,
//                                entity,
//                                entity.hasBoundaryIntersections() ? &(
//                                    *boundary_entity_to_periodic_neighbors_maps_)[type_index][index_set_->index(entity)]
//                                                                  : nullptr);

//  } // ... ibegin(...)

//  IntersectionIterator iend(const typename Codim<0>::Entity& entity) const
//  {
//    const auto& type_index = GlobalGeometryTypeIndex::index(entity.type());
//    assert(!entity.hasBoundaryIntersections()
//           || (*boundary_entity_to_periodic_neighbors_maps_)[type_index].count(index_set_->index(entity)));
//    return IntersectionIterator(BaseType::iend(entity),
//                                *this,
//                                entity,
//                                entity.hasBoundaryIntersections() ? &(
//                                    *boundary_entity_to_periodic_neighbors_maps_)[type_index][index_set_->index(entity)]
//                                                                  : nullptr);
//  } // ... iend(...)

private:
    const MacroElementType& inside_element_;
    const MacroElementType& outside_element_;
    GridGlueType& dd_grid_;
    const MacroGridViewType& macro_grid_view_;
    const std::shared_ptr<LocalGridProviderType> local_inside_grid_;
    std::vector<LocalElementType> inside_elements_;
    std::vector<int> inside_elements_ids_;
    std::vector<std::set<CouplingIntersectionType, CompareType<CouplingIntersectionType>>> coupling_intersections_;
    std::vector<std::pair<size_t, size_t>> local_to_inside_indices_;
}; // ... class CouplingGridViewWrapper ...


} // namespace internal



template <class BaseGridViewImp, class GridGlueType>
class CouplingGridView
  : XT::Common::StorageProvider<internal::CouplingGridViewWrapper<BaseGridViewImp, GridGlueType>>
  , public Dune::GridView<internal::CouplingGridViewWrapperTraits<BaseGridViewImp, GridGlueType>>
{
  static_assert(is_view<BaseGridViewImp>::value);
  using Implementation = internal::CouplingGridViewWrapper<BaseGridViewImp, GridGlueType>;
  using ImplementationStorage = typename XT::Common::StorageProvider<Implementation>;
  using BaseType = Dune::GridView<internal::CouplingGridViewWrapperTraits<BaseGridViewImp, GridGlueType>>;

public:
  using BaseType::dimension;
  using BaseGridViewType = BaseGridViewImp;
  using MacroGridType = typename GridGlueType::MacroGridType;
  using MacroElementType = typename MacroGridType::template Codim<0>::Entity;

  CouplingGridView(const BaseGridViewType& base_grid_view,
                   const MacroElementType& ss,
                   const MacroElementType& nn,
                   GridGlueType& dd_grid)
    : ImplementationStorage(new Implementation(base_grid_view, ss, nn, dd_grid))
    , BaseType(ImplementationStorage::access())
  {}

  CouplingGridView(const CouplingGridView& other)
    : ImplementationStorage(new Implementation(other.access()))
    , BaseType(ImplementationStorage::access())
  {}

  void update()
  {
    ImplementationStorage::access().update();
  }
}; // class CouplingGridView


//template <bool codim_iters_provided, class GL>
//CouplingGridView<GL, codim_iters_provided>
//make_periodic_grid_view(const GL& base_grid_view,
//                        const std::bitset<GL::dimension> periodic_directions = std::bitset<GL::dimension>{}.set())
//{
//  return CouplingGridView<GL, codim_iters_provided>(base_grid_view, periodic_directions);
//}

template <class MG, class GT, class E>
CouplingGridView<MG, GT>
make_coupling_grid_view(const MG& base_grid_view, const E& ss, const E& nn, GT& dd_grid)
{
  return CouplingGridView<MG, GT>(base_grid_view, ss, nn, dd_grid);
}



} // namespace Dune::XT::Grid

#endif // DUNE_XT_GRID_VIEW_PERIODIC_HH
