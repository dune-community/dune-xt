// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Tobias Leibner

#ifndef DUNE_STUFF_GRID_PERIODICVIEW_HH
#define DUNE_STUFF_GRID_PERIODICVIEW_HH

#include <bitset>
#include <map>
#include <utility>
#include <vector>

#if HAVE_DUNE_GRID
#include <dune/grid/common/gridview.hh>
#endif

#include <dune/stuff/aliases.hh>
#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/stuff/common/float_cmp.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/grid/search.hh>

namespace Dune {
namespace Stuff {
namespace Grid {

#if HAVE_DUNE_GRID

namespace internal {

/** \brief Intersection for PeriodicGridView
 *
 * PeriodicIntersection is derived from the Intersection of the underlying GridView. On the inside of the grid or if
 * periodic_ is false, the PeriodicIntersection will behave exactly like its BaseType. If periodic_ is true, the
 * PeriodicIntersection will return neighbor == true even if it actually is on the boundary. In this case, outside(),
 * geometryInOutside() and indexInOutside() are well-defined and give the information from the periodically adjacent
 * entity.
 *
 * \see PeriodicGridView
 */
template <class RealGridViewImp>
class PeriodicIntersection : public RealGridViewImp::Intersection
{
  typedef RealGridViewImp RealGridViewType;
  typedef PeriodicIntersection<RealGridViewType> ThisType;
  typedef typename RealGridViewType::Intersection BaseType;

public:
  using typename BaseType::LocalGeometry;
  typedef typename BaseType::Entity EntityType;
  typedef typename RealGridViewType::IntersectionIterator RealIntersectionIteratorType;
  static const size_t dimDomain = RealGridViewType::dimension;

  //! \brief Constructor from real intersection
  PeriodicIntersection(const BaseType& real_intersection, const RealGridViewType& real_grid_view,
                       const std::pair<bool, EntityType>& periodic_pair)
    : BaseType(real_intersection)
    , periodic_(periodic_pair.first)
    , outside_(periodic_pair.second)
    , real_grid_view_(new DSC::ConstStorageProvider<RealGridViewType>(real_grid_view))
  {
  }

  // methods that differ from BaseType
  bool neighbor() const
  {
    if (periodic_)
      return true;
    else
      return BaseType::neighbor();
  } // bool neighbor() const

  EntityType outside() const
  {
    if (periodic_)
      return outside_;
    else
      return BaseType::outside();
  } // ... outside() const

  LocalGeometry geometryInOutside() const
  {
    if (periodic_) {
      return find_intersection_in_outside().geometryInInside();
    } else {
      return BaseType::geometryInOutside();
    }
  } // ... geometryInOutside() const

  int indexInOutside() const
  {
    if (periodic_) {
      return find_intersection_in_outside().indexInInside();
    } else {
      return BaseType::indexInOutside();
    }
  } // int indexInOutside() const

private:
  // tries to find intersection in outside (works only if periodic_ == true)
  BaseType find_intersection_in_outside() const
  {
    const auto coords                                   = this->geometry().center();
    RealIntersectionIteratorType outside_i_it           = real_grid_view_->access().ibegin(outside_);
    const RealIntersectionIteratorType outside_i_it_end = real_grid_view_->access().iend(outside_);
    // walk over outside intersections and find an intersection on the boundary that differs only in one coordinate
    for (; outside_i_it != outside_i_it_end; ++outside_i_it) {
      const BaseType& curr_outside_intersection = *outside_i_it;
      if (curr_outside_intersection.boundary()) {
        const auto curr_outside_intersection_coords = curr_outside_intersection.geometry().center();
        size_t coord_difference_count = 0;
        for (size_t ii = 0; ii < dimDomain; ++ii) {
          if (Dune::Stuff::Common::FloatCmp::ne(curr_outside_intersection_coords[ii], coords[ii])) {
            ++coord_difference_count;
          }
        }
        if (coord_difference_count == size_t(1)) {
          return *outside_i_it;
        }
      }
    }
    DUNE_THROW(Dune::InvalidStateException, "Could not find outside intersection!");
    return *(real_grid_view_->access().ibegin(outside_));
  } // ... find_intersection_in_outside() const

protected:
  bool periodic_;
  EntityType outside_;
  std::unique_ptr<DSC::ConstStorageProvider<RealGridViewType>> real_grid_view_;
}; // ... class PeriodicIntersection ...

/** \brief IntersectionIterator for PeriodicGridView
 *
 * PeriodicIntersectionIterator is derived from the IntersectionIterator of the underlying GridView and behaves exactly
 * like the underlying IntersectionIterator except that it returns a PeriodicIntersection in its operator* and
 * operator-> methods.
 *
 * \see PeriodicGridView
 */
template <class RealGridViewImp>
class PeriodicIntersectionIterator : public RealGridViewImp::IntersectionIterator
{
  typedef RealGridViewImp RealGridViewType;
  typedef typename RealGridViewType::IntersectionIterator BaseType;

public:
  typedef typename BaseType::Intersection RealIntersectionType;
  typedef int IntersectionIndexType;
  typedef PeriodicIntersection<RealGridViewType> Intersection;
  typedef typename RealGridViewType::template Codim<0>::Entity EntityType;
  typedef std::pair<bool, EntityType> PeriodicPairType;
  static const size_t dimDomain = RealGridViewType::dimension;

  PeriodicIntersectionIterator(BaseType real_intersection_iterator, const RealGridViewType& real_grid_view,
                               const EntityType& entity,
                               const std::map<IntersectionIndexType, PeriodicPairType>& intersection_map)
    : BaseType(real_intersection_iterator)
    , real_grid_view_(real_grid_view)
    , entity_(entity)
    , has_boundary_intersections_(entity_.hasBoundaryIntersections())
    , intersection_map_(intersection_map)
    , nonperiodic_pair_(std::make_pair(bool(false), EntityType(entity_)))
    , current_intersection_(create_current_intersection_safely())
  {
  }

  // methods that differ from BaseType
  const Intersection& operator*() const
  {
    current_intersection_ = create_current_intersection();
    return *current_intersection_;
  }

  const Intersection* operator->() const
  {
    current_intersection_ = create_current_intersection();
    return &(*current_intersection_);
  }

private:
  std::unique_ptr<Intersection> create_current_intersection() const
  {
    return DSC::make_unique<Intersection>(BaseType::operator*(),
                                          real_grid_view_,
                                          has_boundary_intersections_
                                              ? intersection_map_.at((BaseType::operator*()).indexInInside())
                                              : (const PeriodicPairType&)nonperiodic_pair_);
  } // ... create_current_intersection() const

  std::unique_ptr<Intersection> create_current_intersection_safely() const
  {
    const bool is_iend                            = (*this == real_grid_view_.iend(entity_));
    const RealIntersectionType& real_intersection = is_iend ? *real_grid_view_.ibegin(entity_) : BaseType::operator*();
    return DSC::make_unique<Intersection>(real_intersection,
                                          real_grid_view_,
                                          has_boundary_intersections_
                                              ? intersection_map_.at(real_intersection.indexInInside())
                                              : (const PeriodicPairType&)nonperiodic_pair_);
  } // ... create_current_intersection_safely() const

  const RealGridViewType& real_grid_view_;
  const EntityType& entity_;
  const bool has_boundary_intersections_;
  const std::map<IntersectionIndexType, PeriodicPairType>& intersection_map_;
  PeriodicPairType nonperiodic_pair_;
  mutable std::unique_ptr<Intersection> current_intersection_;
}; // ... class PeriodicIntersectionIterator ...

// forward
template <class RealGridViewImp>
class PeriodicGridViewImp;

//! Traits for PeriodicGridView
template <class RealGridViewImp>
class PeriodicGridViewTraits
{
public:
  typedef RealGridViewImp RealGridViewType;
  // use types from RealGridViewType...
  typedef PeriodicGridViewImp<RealGridViewType> GridViewImp;
  typedef typename RealGridViewType::Grid Grid;
  typedef typename RealGridViewType::IndexSet IndexSet;
  typedef typename RealGridViewType::CollectiveCommunication CollectiveCommunication;
  typedef typename RealGridViewType::Traits RealGridViewTraits;

  template <int cd>
  struct Codim
  {
    typedef typename RealGridViewTraits::template Codim<cd>::Iterator Iterator;
    typedef typename RealGridViewTraits::template Codim<cd>::EntityPointer EntityPointer;
    typedef typename RealGridViewTraits::template Codim<cd>::Entity Entity;
    typedef typename RealGridViewTraits::template Codim<cd>::Geometry Geometry;
    typedef typename RealGridViewTraits::template Codim<cd>::LocalGeometry LocalGeometry;

    template <PartitionIteratorType pit>
    struct Partition
    {
      typedef typename RealGridViewTraits::template Codim<cd>::template Partition<pit>::Iterator Iterator;
    };
  }; // ... struct Codim ...

  enum
  {
    conforming = RealGridViewTraits::conforming
  };
  enum
  {
    dimension = Grid::dimension
  };
  enum
  {
    dimensionworld = Grid::dimensionworld
  };

  typedef typename Grid::ctype ctype;

  // ...except for the Intersection and IntersectionIterator
  typedef PeriodicIntersection<RealGridViewType> Intersection;
  typedef PeriodicIntersectionIterator<RealGridViewType> IntersectionIterator;
}; // ... class PeriodicGridViewTraits ...

/** \brief Actual Implementation of PeriodicGridView
 *  \see PeriodicGridView
*/
template <class RealGridViewImp>
class PeriodicGridViewImp : public RealGridViewImp
{
  typedef RealGridViewImp BaseType;
  typedef PeriodicGridViewTraits<BaseType> Traits;

public:
  typedef typename BaseType::Grid Grid;
  typedef typename BaseType::IndexSet IndexSet;
  typedef typename BaseType::template Codim<0>::Entity EntityType;
  typedef PeriodicIntersectionIterator<BaseType> IntersectionIterator;
  typedef typename IntersectionIterator::RealIntersectionType RealIntersectionType;
  typedef typename BaseType::IndexSet::IndexType EntityIndexType;
  typedef int IntersectionIndexType;
  typedef typename RealIntersectionType::GlobalCoordinate DomainType;
  typedef PeriodicIntersection<BaseType> Intersection;
  typedef std::map<IntersectionIndexType, std::pair<bool, EntityType>> IntersectionMapType;
  static const size_t dimDomain = BaseType::dimension;

  template <int cd>
  struct Codim : public Traits::template Codim<cd>
  {
  };

  PeriodicGridViewImp(const BaseType& real_grid_view, const std::bitset<dimDomain> periodic_directions)
    : BaseType(real_grid_view)
    , empty_intersection_map_(IntersectionMapType())
    , periodic_directions_(periodic_directions)
  {
    // find lower left and upper right corner of the grid
    auto entity_it         = BaseType::template begin<0>();
    DomainType lower_left  = entity_it->geometry().center();
    DomainType upper_right = lower_left;
    for (const auto& entity : DSC::entityRange(*this)) {
      const auto i_it_end = BaseType::iend(entity);
      for (auto i_it = BaseType::ibegin(entity); i_it != i_it_end; ++i_it) {
        const RealIntersectionType& intersection = *i_it;
        const auto intersection_coords = intersection.geometry().center();
        for (std::size_t ii = 0; ii < dimDomain; ++ii) {
          if (intersection_coords[ii] > upper_right[ii])
            upper_right[ii] = intersection_coords[ii];
          if (intersection_coords[ii] < lower_left[ii])
            lower_left[ii] = intersection_coords[ii];
        }
      }
    }

    // walk the grid and create a map that maps each entity(index) on the boundary to a map that maps every
    // local intersection index to a std::pair< bool, EntityType >, where the bool component of the pair
    // indicates whether the intersection is on a periodic boundary and the second component is the periodic neighbor
    // entity if the first component is true. If the first component is false, the Entity is not meant to be used.
    EntityInlevelSearch<BaseType> entity_search(*this);
    DomainType periodic_neighbor_coords;
    IntersectionMapType intersection_neighbor_map;
    for (const auto& entity : DSC::entityRange(*this)) {
      if (entity.hasBoundaryIntersections()) {
        intersection_neighbor_map.clear();
        const auto i_it_end = BaseType::iend(entity);
        for (auto i_it = BaseType::ibegin(entity); i_it != i_it_end; ++i_it) {
          const RealIntersectionType& intersection    = *i_it;
          const IntersectionIndexType index_in_inside = intersection.indexInInside();
          bool is_periodic = false;
          if (intersection.boundary()) {
            periodic_neighbor_coords   = intersection.geometry().center();
            size_t num_boundary_coords = 0;
            for (std::size_t ii = 0; ii < dimDomain; ++ii) {
              if (periodic_directions_[ii]) {
                if (Dune::Stuff::Common::FloatCmp::eq(periodic_neighbor_coords[ii], lower_left[ii])) {
                  is_periodic = true;
                  periodic_neighbor_coords[ii] =
                      upper_right[ii] - 1.0 / 100.0 * (entity.geometry().center()[ii] - lower_left[ii]);
                  ++num_boundary_coords;
                } else if (Dune::Stuff::Common::FloatCmp::eq(periodic_neighbor_coords[ii], upper_right[ii])) {
                  is_periodic = true;
                  periodic_neighbor_coords[ii] =
                      lower_left[ii] + 1.0 / 100.0 * (upper_right[ii] - entity.geometry().center()[ii]);
                  ++num_boundary_coords;
                }
              }
            }
            assert(num_boundary_coords = 1);
            if (is_periodic) {
              const auto outside_ptrs = entity_search(std::vector<DomainType>(1, periodic_neighbor_coords));
              const auto& outside_ptr = outside_ptrs.at(0);
              if (outside_ptr == nullptr)
                DUNE_THROW(Dune::InvalidStateException, "Could not find periodic neighbor entity");
              intersection_neighbor_map.insert(
                  std::make_pair(index_in_inside, std::make_pair(is_periodic, *outside_ptr)));
            } else {
              intersection_neighbor_map.insert(
                  std::make_pair(index_in_inside, std::make_pair(is_periodic, EntityType(entity))));
            }
          } else {
            intersection_neighbor_map.insert(
                std::make_pair(index_in_inside, std::make_pair(bool(false), EntityType(entity))));
          }
        }
        entity_to_intersection_map_map_.insert(
            std::make_pair(this->indexSet().index(entity), intersection_neighbor_map));
      }
    }
  } // constructor PeriodicGridViewImp(...)

  IntersectionIterator ibegin(const typename Codim<0>::Entity& entity) const
  {
    return IntersectionIterator(BaseType::ibegin(entity),
                                *this,
                                entity,
                                entity.hasBoundaryIntersections()
                                    ? entity_to_intersection_map_map_.at(this->indexSet().index(entity))
                                    : (const IntersectionMapType&)empty_intersection_map_);

  } // ... ibegin(...)

  IntersectionIterator iend(const typename Codim<0>::Entity& entity) const
  {
    return IntersectionIterator(BaseType::iend(entity),
                                *this,
                                entity,
                                entity.hasBoundaryIntersections()
                                    ? entity_to_intersection_map_map_.at(this->indexSet().index(entity))
                                    : (const IntersectionMapType&)empty_intersection_map_);
  } // ... iend(...)

private:
  std::map<EntityIndexType, IntersectionMapType> entity_to_intersection_map_map_;
  const IntersectionMapType empty_intersection_map_;
  const std::bitset<dimDomain> periodic_directions_;
}; // ... class PeriodicGridViewImp ...

} // namespace internal

/** \brief GridView that takes an arbitrary Dune::GridView and adds periodic boundaries
 *
 * PeriodicGridView is templated by and derived from an arbitrary Dune::GridView. All methods are forwarded to the
 * underlying GridView except for the ibegin and iend methods. These methods return a PeriodicIntersectionIterator
 * which again behaves like the underlying IntersectionIterator except that it returns a PeriodicIntersection in its
 * operator*. The PeriodicIntersection again behaves like an Intersection of the underlying GridView, but may return
 * neighbor() == true and an outside() entity even if it is on the boundary. The outside() entity is the entity
 * adjacent to the intersection if it is identified with the intersection on the other side of the grid.
 * In the constructor, PeriodicGridViewImp will build a map mapping boundary entity indices to a map mapping local
 * intersection indices to a std::pair containing the information whether this intersection shall be periodic and the
 * outside entity. This may take quite long as finding the outside entity requires a grid walk for each periodic
 * intersection.
 * By default, all coordinate directions will be made periodic. By supplying a std::bitset< dimension > you can decide
 * for each direction whether it should be periodic (1 means periodic, 0 means 'behave like underlying GridView in that
 * direction').

   \note
      -  Currently, PeriodicGridView will only work with GridViews on the unit hypercube
      -  Only cube and regular simplex grids have been tested so far. Other grids may not work properly. This is due to
      the heuristics for finding the periodic neighbor entity: Given an intersection on the boundary that shall be
      periodic, the coordinates intersection.geometry().center() are moved to the other side of the grid and then
      supplied to Dune::Stuff::Grid::EntityInLevelSearch. As the coordinates are on the boundary of the wanted entity,
      this search will fail for some grids. Thus, the coordinates are moved a little to the inside of the grid before
      searching for the entity. The moved coordinates will be inside the wanted entity for cube and usual simplex grids
      but this is not guaranteed for arbitrary grids.
 */
template <class RealGridViewImp>
class PeriodicGridView : Dune::Stuff::Common::ConstStorageProvider<internal::PeriodicGridViewImp<RealGridViewImp>>,
                         public Dune::GridView<internal::PeriodicGridViewTraits<RealGridViewImp>>
{
  typedef RealGridViewImp RealGridViewType;
  typedef typename Dune::GridView<internal::PeriodicGridViewTraits<RealGridViewType>> BaseType;
  typedef
      typename Dune::Stuff::Common::ConstStorageProvider<internal::PeriodicGridViewImp<RealGridViewImp>> ConstStorProv;
  typedef typename RealGridViewType::template Codim<0>::Geometry::GlobalCoordinate DomainType;

public:
  static const size_t dimension = RealGridViewType::dimension;

  PeriodicGridView(const RealGridViewType& real_grid_view,
                   const std::bitset<dimension> periodic_directions = std::bitset<dimension>().set())
    : ConstStorProv(new internal::PeriodicGridViewImp<RealGridViewType>(real_grid_view, periodic_directions))
    , BaseType(ConstStorProv::access())
  {
  }

  PeriodicGridView(const PeriodicGridView& other)
    : ConstStorProv(new internal::PeriodicGridViewImp<RealGridViewType>(other.access()))
    , BaseType(ConstStorProv::access())
  {
  }
}; // class PeriodicGridView

#else // HAVE_DUNE_GRID

template <class RealGridViewImp>
class PeriodicGridView
{
  static_assert(AlwaysFalse<RealGridViewImp>::value, "You are missing dune-grid!");
};

#endif // HAVE_DUNE_GRID

} // namespace Grid
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_GRID_PERIODICVIEW_HH
