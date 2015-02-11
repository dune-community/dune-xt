// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Tobias Leibner

#ifndef DUNE_STUFF_GRID_PERIODICVIEW_HH
#define DUNE_STUFF_GRID_PERIODICVIEW_HH

#include <dune/grid/common/gridview.hh>
//#include <dune/grid/common/intersection.hh>
//#include <dune/grid/common/intersectioniterator.hh>

// includes from dune/grid/common/gridview.hh
//#include <dune/common/iteratorrange.hh>
//#include <dune/geometry/type.hh>
//#include <dune/grid/common/datahandleif.hh>
//#include <dune/grid/common/rangegenerators.hh>
// ________________________


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
  /** Copy Constructor from IntersectionImp */
  PeriodicIntersection(const BaseType& real_intersection)
    : BaseType(real_intersection)
  {
  }

  //    //! Copy constructor from an existing intersection.
  //    PeriodicIntersection( PeriodicIntersection& other)
  //      : BaseType(other.real_intersection_)
  //    {}

  //    //! Move constructor from an existing intersection.
  //    PeriodicIntersection(PeriodicIntersection&& other)
  //      : BaseType(std::move(other.real_intersection_))
  //    {}

  /** \brief assignment operator */
  const ThisType& operator=(const BaseType& base)
  {
    real = base.real;
    return *this;
  }


  //! @brief return true if intersection is shared with another element.
  bool neighbor() const
  {
    return real.neighbor();
  }

  using BaseType::inside;

  /*! @brief return EntityPointer to the Entity on the outside of this
     intersection. That is the neighboring Entity.

     @warning Don't call this method if there is no neighboring Entity
     (neighbor() returns false). In this case the result is undefined.
   */
  EntityPointer outside() const
  {
    if (neighbor())
      return ((BaseType*)this)->outside();
    else
      // do fancy stuff
      return ((BaseType*)this)->outside();
  }

  using BaseType::geometryInInside;

  /** \brief geometrical information about this intersection in local
   *         coordinates of the outside() entity.
   *
   *  This method returns a Geometry object that provides a mapping from
   *  local coordinates of the intersection to local coordinates of the
   *  outside() entity.
   *
   *  \note Previously, the geometry was encapsulated in the intersection object
   *        and a const reference was returned.
   *
   *  \note The returned geometry object is guaranteed to remain valid until the
   *        grid is modified (or deleted).
   */
  LocalGeometry geometryInOutside() const
  {
    if (((BaseType*)this)->neighbor())
      return ((BaseType*)this)->geometryInOutside();
    else
      return ((BaseType*)this)->geometryInOutside(); // find periodic neighbor entity and corresponding intersection and
    // use geometryinInside for that Intersection
  }

  using BaseType::indexInInside;

  /** \brief Local index of codim 1 entity in outside() entity where
   *         intersection is contained in
   *
   *  \note This method returns the face number with respect to the generic
   *        reference element.
   *
   *  \returns the index of the outside entity's face containing this
   *           intersection (with respect to the generic reference element)
   */
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
};

template <class RealGridViewImp>
class PeriodicIntersectionIterator : public RealGridViewImp::IntersectionIterator
{

  // type of underlying implementation, for internal use only
  typedef RealGridViewImp GridViewType;
  typedef PeriodicIntersectionIterator<GridViewType> ThisType;
  typedef typename GridViewType::IntersectionIterator BaseType;
  typedef PeriodicIntersection<GridViewType> Intersection;

public:
  /** \brief Type of Intersection this IntersectionIterator points to */

  /** Copy Constructor from IntersectionIteratorImp */
  PeriodicIntersectionIterator(BaseType& real_intersection_iterator)
    : BaseType(real_intersection_iterator)
    , real_intersection_iterator_(real_intersection_iterator)
    , current_intersection_(((BaseType*)this)->operator*())
  {
  }

  PeriodicIntersectionIterator(BaseType real_intersection_iterator)
    : BaseType(real_intersection_iterator)
    , real_intersection_iterator_(real_intersection_iterator)
    , current_intersection_(((BaseType*)this)->operator*())
  {
  }


  //===========================================================
  /** @name Dereferencing
   */
  //@{
  //===========================================================

  /** \brief Dereferencing operator. */
  const Intersection& operator*() const
  {
    return current_intersection_;
  }

  /** \brief Pointer operator. */
  const Intersection* operator->() const
  {
    return new Intersection(((BaseType*)this)->operator*());
  }

  /** @brief Preincrement operator. Proceed to next intersection.*/
  ThisType& operator++()
  {
    BaseType::operator++();
    current_intersection_ = Intersection(BaseType::operator*());
    return *this;
  }

  //===========================================================
  /** @name Compare methods
   */
  //@{
  //===========================================================

  /** @brief Checks for equality.
      Only Iterators pointing to the same intersection from the same Entity
      are equal. Pointing to the same intersection from neighbor is
      unequal as inside and outside are permuted.
   */
  bool operator==(const PeriodicIntersectionIterator& rhs) const
  {
    return rhs.equals(*this);
  }

  /** @brief Checks for inequality.
      Only Iterators pointing to the same intersection from the same Entity
      are equal. Pointing to the same intersection from neighbor is
      unequal as inside and outside are permuted.
   */
  bool operator!=(const PeriodicIntersectionIterator& rhs) const
  {
    return !rhs.equals(*this);
  }
  //@}

  /** @brief forward equality check to realIterator */
  bool equals(const PeriodicIntersectionIterator& rhs) const
  {
    return ((BaseType*)this)->equals(BaseType(rhs));
  }

private:
  Intersection current_intersection_;
  BaseType real_intersection_iterator_;
};

// forward
template <class RealGridViewImp>
class PeriodicGridViewImp;

template <class RealGridViewImp>
class PeriodicGridViewTraits
{
public:
  typedef RealGridViewImp RealGridViewType;
  typedef PeriodicGridViewImp<RealGridViewType> GridViewImp;
  typedef typename RealGridViewType::Grid Grid;
  typedef typename RealGridViewType::IndexSet IndexSet;
  typedef typename RealGridViewType::CollectiveCommunication CollectiveCommunication;

  typedef PeriodicIntersection<RealGridViewType> Intersection;
  typedef PeriodicIntersectionIterator<RealGridViewType> IntersectionIterator;

  typedef typename RealGridViewType::Traits RealGridViewTraits;

  /** \brief A struct that collects all associated types of one implementation
             from the Traits class.
   */
  template <int cd>
  struct Codim
  {
    /** \brief type of iterator returned by the grid view */ // replace this at least for cd == 1 later
    typedef typename RealGridViewTraits::template Codim<cd>::Iterator Iterator;

    /** \brief type of corresponding entity pointer */
    typedef typename RealGridViewTraits::template Codim<cd>::EntityPointer EntityPointer;

    /** \brief type of corresponding entity */
    typedef typename RealGridViewTraits::template Codim<cd>::Entity Entity;

    /** \brief type of the geometry implementation */
    typedef typename RealGridViewTraits::template Codim<cd>::Geometry Geometry;

    /** \brief type of the implementation for local geometries */
    typedef typename RealGridViewTraits::template Codim<cd>::LocalGeometry LocalGeometry;

    /** \brief Define types needed to iterate over entities of a given partition type */
    template <PartitionIteratorType pit>
    struct Partition
    {
      /** \brief iterator over a given codim and partition type */
      typedef typename RealGridViewTraits::template Codim<cd>::template Partition<pit>::Iterator Iterator;
    };
  }; //: public Traits :: template Codim<cd> {};

  enum
  {
    /** \brief Export if this grid view is conforming */
    conforming = RealGridViewTraits::conforming
  };

  /** \brief type used for coordinates in grid */
  typedef typename Grid::ctype ctype;

  enum
  { //! \brief The dimension of the grid
    dimension = Grid::dimension
  };

  enum
  { //! \brief The dimension of the world the grid lives in
    dimensionworld = Grid::dimensionworld
  };
};


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

  /** \brief Codim Structure */
  template <int cd>
  struct Codim : public Traits::template Codim<cd>
  {
  };

  /** \brief constructor (engine concept) */
  PeriodicGridViewImp(const RealGridViewType& real_grid_view)
    : real_grid_view_(real_grid_view)
  {
  }

  /** \brief obtain a const reference to the underlying hierarchic grid */
  const Grid& grid() const
  {
    return real_grid_view_.grid();
  }

  /** \brief obtain the index set
   *
   * The lifetime of the returned index set is bound to the lifetime of the
   * grid view. Keep a copy of the grid view to prevent the index set from
   * becoming a dangling reference.
   */
  const IndexSet& indexSet() const
  {
    return real_grid_view_.indexSet();
  }

  /** \brief obtain number of entities in a given codimension */
  int size(int codim) const
  {
    return real_grid_view_.size(codim);
  }

  /** \brief obtain number of entities with a given geometry type */
  int size(const GeometryType& type) const
  {
    return real_grid_view_.size(type);
  }

  /** @brief Return true if the given entity is contained in this grid view
   * @todo Currently we call the implementation on the IndexSet.  This may lead to suboptimal efficiency.
   *
   * \note If the input element e is not an element of the grid, then
   *       the result of contains() is undefined.
   */
  template <class EntityType>
  bool contains(const EntityType& e) const
  {
    return real_grid_view_.indexSet().contains(e);
  }

  /** \brief obtain begin iterator for this view */
  template <int cd>
  typename Codim<cd>::Iterator begin() const
  {
    return real_grid_view_.template begin<cd>();
  }

  /** \brief obtain end iterator for this view */
  template <int cd>
  typename Codim<cd>::Iterator end() const
  {
    return real_grid_view_.template end<cd>();
  }

  /** \brief obtain begin iterator for this view */
  template <int cd, PartitionIteratorType pitype>
  typename Codim<cd>::template Partition<pitype>::Iterator begin() const
  {
    return real_grid_view_.template begin<cd, pitype>();
  }

  /** \brief obtain end iterator for this view */
  template <int cd, PartitionIteratorType pitype>
  typename Codim<cd>::template Partition<pitype>::Iterator end() const
  {
    return real_grid_view_.template end<cd, pitype>();
  }

  /** \brief obtain begin intersection iterator with respect to this view */
  IntersectionIterator ibegin(const typename Codim<0>::Entity& entity) const
  {
    return IntersectionIterator(real_grid_view_.ibegin(entity));
  }

  /** \brief obtain end intersection iterator with respect to this view */
  IntersectionIterator iend(const typename Codim<0>::Entity& entity) const
  {
    return IntersectionIterator(real_grid_view_.iend(entity));
  }

  /** \brief obtain collective communication object */
  const CollectiveCommunication& comm() const
  {
    return real_grid_view_.comm();
  }

  /** \brief Return size of the overlap region for a given codim on the grid view.  */
  int overlapSize(int codim) const
  {
    return real_grid_view_.overlapSize(codim);
  }

  /** \brief Return size of the ghost region for a given codim on the grid view.  */
  int ghostSize(int codim) const
  {
    return real_grid_view_.ghostSize(codim);
  }

  /** \brief Communicate data on this view */
  template <class DataHandleImp, class DataType>
  void communicate(CommDataHandleIF<DataHandleImp, DataType>& data, InterfaceType iftype,
                   CommunicationDirection dir) const
  {
    real_grid_view_.communicate(data, iftype, dir);
  }

private:
  const RealGridViewType& real_grid_view_;
};

template <class RealGridViewImp>
class PeriodicGridView : public Dune::GridView<PeriodicGridViewTraits<RealGridViewImp>>
{
  typedef RealGridViewImp RealGridViewType;
  typedef typename Dune::GridView<PeriodicGridViewTraits<RealGridViewType>> BaseType;
  typedef PeriodicGridViewImp<RealGridViewType> PeriodicDings;

public:
  PeriodicGridView(const RealGridViewType& real_grid_view)
    : BaseType(*(new PeriodicDings(real_grid_view)))
    , periodic_grid_view_(real_grid_view)
  {
  }

private:
  const PeriodicGridViewImp<RealGridViewType>& periodic_grid_view_;
};

} // namespace Grid
} // namespace Stuff
} // namespace Dune

#endif // #ifndef DUNE_STUFF_GRID_PERIODICVIEW_HH
