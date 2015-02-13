// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Tobias Leibner

#ifndef DUNE_STUFF_GRID_PERIODICVIEW_HH
#define DUNE_STUFF_GRID_PERIODICVIEW_HH

#include <bitset>
#include <vector>

#include <dune/grid/common/gridview.hh>

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
                       const std::bitset<dimDomain> periodic_directions)
    : BaseType(real_intersection)
    , periodic_(false)
    , real_grid_view_(&real_grid_view)
    , periodic_neighbor_coords_(1)
    , num_boundary_coords_(0)
  {
    if (BaseType::boundary()) {
      periodic_neighbor_coords_[0] = BaseType::geometry().center();
      for (std::size_t ii = 0; ii < dimDomain; ++ii) {
        if (Dune::Stuff::Common::FloatCmp::eq(periodic_neighbor_coords_[0][ii], 0.0)) {
          if (periodic_directions[ii]) {
            periodic_                        = true;
            periodic_neighbor_coords_[0][ii] = 1.0;
            ++num_boundary_coords_;
          }
        } else if (Dune::Stuff::Common::FloatCmp::eq(periodic_neighbor_coords_[0][ii], 1.0)) {
          if (periodic_directions[ii]) {
            periodic_                        = true;
            periodic_neighbor_coords_[0][ii] = 0.0;
            ++num_boundary_coords_;
          }
        }
      }
      assert(num_boundary_coords_ = 1);
    }
  } // constructor

  //! \brief Copy operator from existing PeriodicIntersection
  const ThisType operator=(ThisType other)
  {
    real                      = other.real;
    periodic_                 = other.periodic_;
    periodic_neighbor_coords_ = other.periodic_neighbor_coords_;
    real_grid_view_           = other.real_grid_view_;
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
      return *(Dune::Stuff::Grid::EntityInlevelSearch<RealGridViewType>(*real_grid_view_)
                   .
                   operator()(periodic_neighbor_coords_)[0]);
    else
      return BaseType::outside();
  }

  LocalGeometry geometryInOutside() const
  {
    if (periodic_) {
      auto& outside_entity = outside();
      for (auto i_it = real_grid_view_->ibegin(outside_entity); i_it != real_grid_view_->iend(outside_entity); ++i_it) {
        auto& intersection = *i_it;
        if (Dune::Stuff::Common::FloatCmp::eq(intersection.geometry().center(), periodic_neighbor_coords_[0]))
          return intersection.geometryInInside();
      }
    } else
      return BaseType::geometryInOutside();
  } // geometryInOutSide

  int indexInOutside() const
  {
    if (periodic_) {
      auto& outside_entity = outside();
      for (auto i_it = *real_grid_view_->ibegin(outside_entity); i_it != real_grid_view_->iend(outside_entity);
           ++i_it) {
        auto& intersection = *i_it;
        if (Dune::Stuff::Common::FloatCmp::eq(intersection.geometry().center(), periodic_neighbor_coords_[0]))
          return intersection.indexInInside();
      }
    } else
      return BaseType::indexInOutside();
  }

protected:
  using BaseType::real;
  bool periodic_;
  std::vector<GlobalCoordinate> periodic_neighbor_coords_;
  const RealGridViewType* real_grid_view_;
  std::size_t num_boundary_coords_;
}; // ... class PeriodicIntersection ...

template <class RealGridViewImp>
class PeriodicIntersectionIterator : public RealGridViewImp::IntersectionIterator
{
  typedef RealGridViewImp RealGridViewType;
  typedef PeriodicIntersectionIterator<RealGridViewType> ThisType;
  typedef typename RealGridViewType::IntersectionIterator BaseType;
  typedef PeriodicIntersection<RealGridViewType> Intersection;
  typedef typename RealGridViewType::template Codim<0>::Entity EntityType;
  static const size_t dimDomain = RealGridViewType::dimension;

public:
  /** Copy Constructor from real intersection iterator */
  PeriodicIntersectionIterator(BaseType& real_intersection_iterator, const RealGridViewType& real_grid_view,
                               const std::bitset<dimDomain> periodic_directions, const EntityType& entity)
    : BaseType(real_intersection_iterator)
    , current_intersection_((real_intersection_iterator == real_grid_view.iend(entity)) ? *real_grid_view.ibegin(entity)
                                                                                        : BaseType::
                                                                                          operator*(),
                            real_grid_view, periodic_directions)
    , periodic_directions_(periodic_directions)
    , real_grid_view_(real_grid_view)
    , entity_(entity)
  {
  }

  /** Copy Constructor from real intersection iterator, pass by value */
  PeriodicIntersectionIterator(BaseType real_intersection_iterator, const RealGridViewType& real_grid_view,
                               const std::bitset<dimDomain> periodic_directions, const EntityType& entity)
    : BaseType(real_intersection_iterator)
    , current_intersection_((real_intersection_iterator == real_grid_view.iend(entity)) ? *real_grid_view.ibegin(entity)
                                                                                        : BaseType::
                                                                                          operator*(),
                            real_grid_view, periodic_directions)
    , periodic_directions_(periodic_directions)
    , real_grid_view_(real_grid_view)
    , entity_(entity)
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
        periodic_directions_);
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
    return BaseType::equals(BaseType(rhs)) && (current_intersection_ == rhs.current_intersection_);
  }

private:
  Intersection current_intersection_;
  std::bitset<dimDomain> periodic_directions_;
  const RealGridViewType& real_grid_view_;
  const EntityType& entity_;
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
  typedef PeriodicIntersection<RealGridViewType> Intersection;
  typedef typename RealGridViewType::CollectiveCommunication CollectiveCommunication;
  static const size_t dimDomain = RealGridViewType::dimension;

  template <int cd>
  struct Codim : public Traits::template Codim<cd>
  {
  };

  /** Constructor from real grid view */
  PeriodicGridViewImp(const RealGridViewType& real_grid_view, const std::bitset<dimDomain> periodic_directions)
    : real_grid_view_(real_grid_view)
    , periodic_directions_(periodic_directions)
  {
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
    return IntersectionIterator(real_grid_view_.ibegin(entity), real_grid_view_, periodic_directions_, entity);
  }

  IntersectionIterator iend(const typename Codim<0>::Entity& entity) const
  {
    return IntersectionIterator(real_grid_view_.iend(entity), real_grid_view_, periodic_directions_, entity);
  }

private:
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
  static const size_t dimDomain = RealGridViewType::dimension;

public:
  PeriodicGridView(const RealGridViewType& real_grid_view,
                   const std::bitset<dimDomain> periodic_directions = std::bitset<dimDomain>().set())
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
