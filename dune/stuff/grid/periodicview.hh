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

#include <dune/grid/common/gridview.hh>

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/stuff/common/float_cmp.hh>
#include <dune/stuff/grid/search.hh>

#include <dune/stuff/common/string.hh>

namespace Dune {
namespace Stuff {
namespace Grid {

template <class RealGridViewImp>
class PeriodicIntersection : public RealGridViewImp::Intersection
{
  typedef RealGridViewImp RealGridViewType;
  typedef PeriodicIntersection<RealGridViewType> ThisType;
  typedef typename RealGridViewType::Intersection BaseType;

  typedef typename BaseType::LocalGeometry LocalGeometry;
  typedef typename BaseType::Geometry::GlobalCoordinate GlobalCoordinate;
  typedef typename BaseType::EntityPointer EntityPointer;
  static const size_t dimDomain = RealGridViewType::dimension;

public:
  //! \brief Constructor from real intersection
  PeriodicIntersection(const BaseType& real_intersection, const RealGridViewType& real_grid_view,
                       const std::tuple<bool, EntityPointer, const BaseType*>& periodic_tuple)
    : BaseType(real_intersection)
    , periodic_(std::get<0>(periodic_tuple))
    , outside_(std::get<1>(periodic_tuple))
    , outside_intersection_(std::get<2>(periodic_tuple))
    , real_grid_view_(&real_grid_view)
  {
  }

  //! \brief Copy operator from existing PeriodicIntersection
  const ThisType operator=(ThisType other)
  {
    real                  = other.real;
    periodic_             = other.periodic_;
    outside_              = other.outside_;
    outside_intersection_ = other.outside_intersection_;
    real_grid_view_       = other.real_grid_view_;
    return *this;
  }

  // methods that differ from BaseType
  bool neighbor() const
  {
    if (periodic_)
      return true;
    else
      return BaseType::neighbor();
  }

  bool boundary() const
  {
    if (periodic_) {
      return false;
    } else {
      if (BaseType::boundary()) {
        std::cout << "Komisch " << BaseType::geometry().center() << std::endl;
      }
      return BaseType::boundary();
    }
  }

  EntityPointer outside() const
  {
    if (periodic_)
      return outside_;
    else
      return BaseType::outside();
  }

  LocalGeometry geometryInOutside() const
  {
    if (periodic_)
      return outside_intersection_->geometryInInside();
    else
      return BaseType::geometryInOutside();
  } // geometryInOutSide

  int indexInOutside() const
  {
    if (periodic_)
      return outside_intersection_->indexInInside();
    else
      return BaseType::indexInOutside();
  }

protected:
  using BaseType::real;
  bool periodic_;
  EntityPointer outside_;
  const RealGridViewType* real_grid_view_;
  const BaseType* outside_intersection_;
}; // ... class PeriodicIntersection ...

template <class RealGridViewImp>
class PeriodicIntersectionIterator : public RealGridViewImp::IntersectionIterator
{
  typedef RealGridViewImp RealGridViewType;
  typedef PeriodicIntersectionIterator<RealGridViewType> ThisType;
  typedef typename RealGridViewType::IntersectionIterator BaseType;

public:
  typedef typename BaseType::Intersection RealIntersectionType;
  typedef int IntersectionIndexType;
  typedef typename RealIntersectionType::EntityPointer EntityPointerType;
  typedef PeriodicIntersection<RealGridViewType> Intersection;
  typedef typename RealGridViewType::template Codim<0>::Entity EntityType;
  static const size_t dimDomain = RealGridViewType::dimension;

  /** Copy Constructor from real intersection iterator */
  PeriodicIntersectionIterator(
      BaseType& real_intersection_iterator, const RealGridViewType& real_grid_view, const EntityType& entity,
      const std::map<IntersectionIndexType, std::tuple<bool, EntityPointerType, const RealIntersectionType*>>&
          intersection_map)
    : BaseType(real_intersection_iterator)
    , current_intersection_(
          (real_intersection_iterator == real_grid_view.iend(entity)) ? *real_grid_view.ibegin(entity) : BaseType::
                                                                                                         operator*(),
          real_grid_view, intersection_map.at((real_intersection_iterator == real_grid_view.iend(entity))
                                                  ? *real_grid_view.ibegin(entity)
                                                  : BaseType::
                                                    operator*())
                              .indexInInside)
    , real_grid_view_(real_grid_view)
    , entity_(entity)
    , intersection_map_(intersection_map)
  {
  }

  /** Copy Constructor from real intersection iterator, pass by value */
  PeriodicIntersectionIterator(
      BaseType real_intersection_iterator, const RealGridViewType& real_grid_view, const EntityType& entity,
      const std::map<IntersectionIndexType, std::tuple<bool, EntityPointerType, const RealIntersectionType*>>&
          intersection_map)
    : BaseType(real_intersection_iterator)
    , current_intersection_(
          (real_intersection_iterator == real_grid_view.iend(entity)) ? *real_grid_view.ibegin(entity) : BaseType::
                                                                                                         operator*(),
          real_grid_view, intersection_map.at(((real_intersection_iterator == real_grid_view.iend(entity))
                                                   ? *real_grid_view.ibegin(entity)
                                                   : BaseType::
                                                     operator*())
                                                  .indexInInside()))
    , real_grid_view_(real_grid_view)
    , entity_(entity)
    , intersection_map_(intersection_map)
  {
  }


  // methods that differ from BaseType
  const Intersection& operator*() const
  {
    return current_intersection_;
  }

  const Intersection* operator->() const
  {
    return &current_intersection_;
  }

  ThisType& operator++()
  {
    BaseType::operator++();
    current_intersection_ = Intersection(
        BaseType::equals(real_grid_view_.iend(entity_)) ? *real_grid_view_.ibegin(entity_) : BaseType::operator*(),
        real_grid_view_,
        intersection_map_.at(
            (BaseType::equals(real_grid_view_.iend(entity_)) ? *real_grid_view_.ibegin(entity_) : BaseType::operator*())
                .indexInInside()));
    return *this;
  }

  bool operator==(const PeriodicIntersectionIterator& rhs) const
  {
    return rhs.equals(*this);
  }

  bool operator!=(const PeriodicIntersectionIterator& rhs) const
  {
    return !rhs.equals(*this);
  }

  bool equals(const PeriodicIntersectionIterator& rhs) const
  {
    return BaseType::equals(BaseType(rhs)) /*&& (current_intersection_ == rhs.current_intersection_)*/;
  }

private:
  Intersection current_intersection_;
  const RealGridViewType& real_grid_view_;
  const EntityType& entity_;
  const std::map<IntersectionIndexType, std::tuple<bool, EntityPointerType, const RealIntersectionType*>>&
      intersection_map_;

}; // ... class PeriodicIntersectionIterator ...

// forward
template <class RealGridViewImp>
class PeriodicGridViewImp;

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
    // TODO: replace this at least for cd == 1 later?
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


template <class RealGridViewImp>
class PeriodicGridViewImp
{
  typedef RealGridViewImp RealGridViewType;
  typedef PeriodicGridViewImp<RealGridViewType> ThisType;
  typedef PeriodicGridViewTraits<RealGridViewType> Traits;

public:
  typedef typename RealGridViewType::Grid Grid;
  typedef typename RealGridViewType::IndexSet IndexSet;
  typedef typename Dune::GeometryType GeometryType;
  typedef PeriodicIntersectionIterator<RealGridViewType> IntersectionIterator;
  typedef typename IntersectionIterator::RealIntersectionType RealIntersectionType;
  typedef typename RealGridViewType::IndexSet::IndexType EntityIndexType;
  typedef int IntersectionIndexType;
  typedef typename RealIntersectionType::GlobalCoordinate CoordinateType;
  typedef PeriodicIntersection<RealGridViewType> Intersection;
  typedef typename RealGridViewType::CollectiveCommunication CollectiveCommunication;
  static const size_t dimDomain = RealGridViewType::dimension;

  typedef typename Grid::template Codim<0>::EntityPointer EntityPointerType;

  template <int cd>
  struct Codim : public Traits::template Codim<cd>
  {
  };

  /** Constructor from real grid view */
  PeriodicGridViewImp(const RealGridViewType& real_grid_view, const std::bitset<dimDomain> periodic_directions)
    : real_grid_view_(real_grid_view)
    , periodic_directions_(periodic_directions)
  {
    const auto& it_end = real_grid_view_.template end<0>();
    for (auto it = real_grid_view_.template begin<0>(); it != it_end; ++it) {
      const auto& entity = *it;
      std::map<IntersectionIndexType, std::tuple<bool, EntityPointerType, const RealIntersectionType*>>
          intersection_neighbor_map;
      //      if (entity.partitionType() == BorderEntity) {
      const auto& i_it_end = real_grid_view_.iend(entity);
      for (auto i_it = real_grid_view_.ibegin(entity); i_it != i_it_end; ++i_it) {
        const RealIntersectionType& intersection = *i_it;
        bool is_periodic = false;
        if (intersection.boundary()) {
          CoordinateType periodic_neighbor_coords = intersection.geometry().center();
          size_t num_boundary_coords = 0;
          for (std::size_t ii = 0; ii < dimDomain; ++ii) {
            if (Dune::Stuff::Common::FloatCmp::eq(periodic_neighbor_coords[ii], 0.0)) {
              if (periodic_directions_[ii]) {
                is_periodic                  = true;
                periodic_neighbor_coords[ii] = 1.0;
                ++num_boundary_coords;
              }
            } else if (Dune::Stuff::Common::FloatCmp::eq(periodic_neighbor_coords[ii], 1.0)) {
              if (periodic_directions_[ii]) {
                is_periodic                  = true;
                periodic_neighbor_coords[ii] = 0.0;
                ++num_boundary_coords;
              }
            }
          }
          assert(num_boundary_coords = 1);
          if (is_periodic) {
            EntityPointerType periodic_neighbor =
                *(Dune::Stuff::Grid::EntityInlevelSearch<RealGridViewType>(real_grid_view_)
                      .
                      operator()(std::vector<CoordinateType>(1, periodic_neighbor_coords))[0]);
            const RealIntersectionType* outside_intersection_pointer;
            const auto& neighbor_i_it_end = real_grid_view_.iend(*periodic_neighbor);
            for (auto neighbor_i_it = real_grid_view_.ibegin(*periodic_neighbor); neighbor_i_it != neighbor_i_it_end;
                 ++neighbor_i_it) {
              const auto& outside_intersection = *i_it;
              if (Dune::Stuff::Common::FloatCmp::eq(periodic_neighbor_coords,
                                                    outside_intersection.geometry().center())) {
                outside_intersection_pointer = &outside_intersection;
              }
            }
            const size_t index_in_inside = intersection.indexInInside();
            const auto pair = std::make_pair(
                index_in_inside, std::make_tuple(is_periodic, periodic_neighbor, outside_intersection_pointer));
            intersection_neighbor_map.insert(pair);
          } else {
            intersection_neighbor_map.insert(std::make_pair(
                intersection.indexInInside(), std::make_tuple(is_periodic, EntityPointerType(entity), &intersection)));
          }
        } else {
          intersection_neighbor_map.insert(std::make_pair(
              intersection.indexInInside(), std::make_tuple(bool(false), EntityPointerType(entity), &intersection)));
        }
      }
      entity_to_intersection_map_map_.insert(
          std::make_pair(real_grid_view_.indexSet().index(entity), intersection_neighbor_map));
      //      }
    }
  }

  // forward methods to real grid view ...
  const Grid& grid() const
  {
    return real_grid_view_.grid();
  }

  const IndexSet& indexSet() const
  {
    return real_grid_view_.indexSet();
  }

  int size(int codim) const
  {
    return real_grid_view_.size(codim);
  }

  int size(const GeometryType& type) const
  {
    return real_grid_view_.size(type);
  }

  template <class EntityType>
  bool contains(const EntityType& e) const
  {
    return real_grid_view_.indexSet().contains(e);
  }

  template <int cd>
  typename Codim<cd>::Iterator begin() const // TODO: change at least for cd == 1 to PeriodicIterator
  {
    return real_grid_view_.template begin<cd>();
  }

  template <int cd>
  typename Codim<cd>::Iterator end() const
  {
    return real_grid_view_.template end<cd>();
  }

  template <int cd, PartitionIteratorType pitype>
  typename Codim<cd>::template Partition<pitype>::Iterator begin() const
  {
    return real_grid_view_.template begin<cd, pitype>();
  }

  template <int cd, PartitionIteratorType pitype>
  typename Codim<cd>::template Partition<pitype>::Iterator end() const
  {
    return real_grid_view_.template end<cd, pitype>();
  }

  const CollectiveCommunication& comm() const
  {
    return real_grid_view_.comm();
  }

  int overlapSize(int codim) const
  {
    return real_grid_view_.overlapSize(codim);
  }

  int ghostSize(int codim) const
  {
    return real_grid_view_.ghostSize(codim);
  }

  template <class DataHandleImp, class DataType>
  void communicate(CommDataHandleIF<DataHandleImp, DataType>& data, InterfaceType iftype,
                   CommunicationDirection dir) const
  {
    real_grid_view_.communicate(data, iftype, dir);
  }

  // ... except for the intersection iterators
  IntersectionIterator ibegin(const typename Codim<0>::Entity& entity) const
  {
    return IntersectionIterator(real_grid_view_.ibegin(entity),
                                real_grid_view_,
                                entity,
                                entity_to_intersection_map_map_.at(real_grid_view_.indexSet().index(entity)));
  }

  IntersectionIterator iend(const typename Codim<0>::Entity& entity) const
  {
    return IntersectionIterator(real_grid_view_.iend(entity),
                                real_grid_view_,
                                entity,
                                entity_to_intersection_map_map_.at(real_grid_view_.indexSet().index(entity)));
  }

private:
  std::map<EntityIndexType,
           std::map<IntersectionIndexType, std::tuple<bool, EntityPointerType, const RealIntersectionType*>>>
      entity_to_intersection_map_map_;
  const RealGridViewType& real_grid_view_;
  const std::bitset<dimDomain> periodic_directions_;
}; // ... class PeriodicGridViewImp ...

template <class RealGridViewImp>
class PeriodicGridView : Dune::Stuff::Common::ConstStorageProvider<PeriodicGridViewImp<RealGridViewImp>>,
                         public Dune::GridView<PeriodicGridViewTraits<RealGridViewImp>>
{
  typedef RealGridViewImp RealGridViewType;
  typedef typename Dune::GridView<PeriodicGridViewTraits<RealGridViewType>> BaseType;
  typedef typename Dune::Stuff::Common::ConstStorageProvider<PeriodicGridViewImp<RealGridViewImp>> ConstStorProv;

public:
  static const size_t dimension = RealGridViewType::dimension;

  PeriodicGridView(const RealGridViewType& real_grid_view,
                   const std::bitset<dimension> periodic_directions = std::bitset<dimension>().set())
    : ConstStorProv(new PeriodicGridViewImp<RealGridViewType>(real_grid_view, periodic_directions))
    , BaseType(ConstStorProv::access())
  {
  }

  PeriodicGridView(const PeriodicGridView& other)
    : ConstStorProv(other.access())
    , BaseType(ConstStorProv::access())
  {
  }
}; // ... class PeriodicGridView ...

} // namespace Grid
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_GRID_PERIODICVIEW_HH
