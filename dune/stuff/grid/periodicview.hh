// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Tobias Leibner

#ifndef DUNE_STUFF_GRID_PERIODICVIEW_HH
#define DUNE_STUFF_GRID_PERIODICVIEW_HH

#include <dune/grid/common/gridview.hh>

#include <dune/stuff/common/memory.hh>

namespace Dune {
namespace Stuff {
namespace Grid {

template <class RealGridViewImp>
class PeriodicIntersection : public RealGridViewImp::Intersection
{
  typedef PeriodicIntersection<RealGridViewImp> ThisType;
  typedef typename RealGridViewImp::Intersection BaseType;

  typedef typename BaseType::LocalGeometry LocalGeometry;
  typedef typename BaseType::EntityPointer EntityPointer;

public:
  //! \brief Constructor from real intersection
  PeriodicIntersection(const BaseType& real_intersection)
    : BaseType(real_intersection)
  {
  }

  //! \brief Copy operator from existing PeriodicIntersection
  const ThisType operator=(ThisType other)
  {
    real = other.real;
    return *this;
  }

  // methods that differ from BaseType
  bool neighbor() const
  {
    return real.neighbor();
  }

  EntityPointer outside() const
  {
    if (neighbor())
      return ((BaseType*)this)->outside();
    else
      // do fancy stuff
      return ((BaseType*)this)->outside();
  }

  LocalGeometry geometryInOutside() const
  {
    if (((BaseType*)this)->neighbor())
      return ((BaseType*)this)->geometryInOutside();
    else
      return ((BaseType*)this)->geometryInOutside(); // find periodic neighbor entity and corresponding intersection and
    // use geometryinInside for that Intersection
  }

  int indexInOutside() const
  {
    if (((BaseType*)this)->neighbor())
      return ((BaseType*)this)->indexInOutside();
    else
      return ((BaseType*)this)->indexInOutside(); // find periodic neighbor entity and corresponding intersection and
    // use indexInInside for that Intersection
  }

protected:
  using BaseType::real;
}; // ... class PeriodicIntersection ...

template <class RealGridViewImp>
class PeriodicIntersectionIterator : public RealGridViewImp::IntersectionIterator
{
  typedef RealGridViewImp GridViewType;
  typedef PeriodicIntersectionIterator<GridViewType> ThisType;
  typedef typename GridViewType::IntersectionIterator BaseType;
  typedef PeriodicIntersection<GridViewType> Intersection;

public:
  /** Copy Constructor from real intersection iterator */
  PeriodicIntersectionIterator(BaseType& real_intersection_iterator)
    : BaseType(real_intersection_iterator)
    , current_intersection_(((BaseType*)this)->operator*())
  {
  }

  /** Copy Constructor from real intersection iterator, pass by value */
  PeriodicIntersectionIterator(BaseType real_intersection_iterator)
    : BaseType(real_intersection_iterator)
    , current_intersection_(((BaseType*)this)->operator*())
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
    current_intersection_ = Intersection(BaseType::operator*());
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
    return ((BaseType*)this)->equals(BaseType(rhs));
  }

private:
  Intersection current_intersection_;
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

  template <int cd>
  struct Codim : public Traits::template Codim<cd>
  {
  };

  /** Constructor from real grid view */
  PeriodicGridViewImp(const RealGridViewType& real_grid_view)
    : real_grid_view_(real_grid_view)
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
    return IntersectionIterator(real_grid_view_.ibegin(entity));
  }

  IntersectionIterator iend(const typename Codim<0>::Entity& entity) const
  {
    return IntersectionIterator(real_grid_view_.iend(entity));
  }

private:
  const RealGridViewType& real_grid_view_;
}; // ... class PeriodicGridViewImp ...

template <class RealGridViewImp>
class PeriodicGridView : Dune::Stuff::Common::ConstStorageProvider<PeriodicGridViewImp<RealGridViewImp>>,
                         public Dune::GridView<PeriodicGridViewTraits<RealGridViewImp>>
{
  typedef RealGridViewImp RealGridViewType;
  typedef typename Dune::GridView<PeriodicGridViewTraits<RealGridViewType>> BaseType;
  typedef typename Dune::Stuff::Common::ConstStorageProvider<PeriodicGridViewImp<RealGridViewImp>> ConstStorProv;

public:
  PeriodicGridView(const RealGridViewType& real_grid_view)
    : ConstStorProv(new PeriodicGridViewImp<RealGridViewType>(real_grid_view))
    , BaseType(ConstStorProv::access())
  {
  }

  PeriodicGridView(const PeriodicGridView& other)
    : ConstStorProv(other.access())
    , BaseType(ConstStorProv::access())
  {
  }
};

} // namespace Grid
} // namespace Stuff
} // namespace Dune

#endif // #ifndef DUNE_STUFF_GRID_PERIODICVIEW_HH
