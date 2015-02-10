// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Tobias Leibner

#ifndef DUNE_STUFF_GRID_PERIODICVIEW_HH
#define DUNE_STUFF_GRID_PERIODICVIEW_HH

#include <dune/grid/common/gridview.hh>
#include <dune/grid/common/intersection.hh>
#include <dune/grid/common/intersectioniterator.hh>

// includes from dune/grid/common/gridview.hh
#include <dune/common/iteratorrange.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/rangegenerators.hh>
// ________________________


namespace Dune {
namespace Stuff {
namespace Grid {

template <class GridViewImp>
class PeriodicIntersection : GridViewImp::Intersection
{
  typedef PeriodicIntersection<GridViewImp> ThisType;
  typedef GridViewImp::Intersection BaseType;
  //  typedef typename GridViewImp::Traits::Grid   GridType;
public:
  /** \brief Type of entity that this Intersection belongs to */
  typedef typename BaseType::template Codim<0>::Entity Entity;

  /** \brief Pointer to the type of entities that this Intersection belongs to */
  typedef typename BaseType::template Codim<0>::EntityPointer EntityPointer;

  /** \brief Codim 1 geometry returned by geometry() */
  typedef typename BaseType::template Codim<1>::Geometry Geometry;

  /** \brief Type for vectors of coordinates on the intersection */
  typedef typename Geometry::LocalCoordinate LocalCoordinate;

  /** \brief Type for normal vectors */
  typedef typename Geometry::GlobalCoordinate GlobalCoordinate;

  /** \brief Codim 1 geometry returned by geometryInInside() and geometryInOutside() */
  typedef typename BaseType::template Codim<1>::LocalGeometry LocalGeometry;

  //! Type of individual coefficients of coordinate vectors
  typedef typename BaseType::ctype ctype;

  /** Copy Constructor from IntersectionImp */
  PeriodicIntersection(const BaseType& real_intersection)
    : BaseType(real_intersection)
    , real_intersection_(real_intersection)
  {
  }

  //! Copy constructor from an existing intersection.
  PeriodicIntersection(PeriodicIntersection& other)
    : BaseType(other.real_intersection_)
    , real_intersection_(other.real_intersection_)
  {
  }

  //! Move constructor from an existing intersection.
  PeriodicIntersection(PeriodicIntersection&& other)
    : BaseType(std::move(other.real_intersection_))
    , real_intersection_(std::move(other.real_intersection_))
  {
  }

  //! Copy assignment operator from an existing intersection.
  PeriodicIntersection& operator=(const PeriodicIntersection& other)
  {
    real_intersection_ = other.real_intersection_;
    return *this;
  }

  //! Move assignment operator from an existing intersection.
  PeriodicIntersection& operator=(PeriodicIntersection&& other)
  {
    real_intersection_ = std::move(other.real_intersection_);
    return *this;
  }

  //! Return true if intersection is with interior or exterior boundary (see the cases above)
  bool boundary() const
  {
    return real_intersection_.boundary();
  }

#if DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS
  /**
     \brief Identifier for boundary segment from macro grid.

     One can attach a boundary Id to a boundary segment on the macro
     grid. This Id will also be used for all fragments of these
     boundary segments.

     The numbering is defined as:
     - Id==0 for all intersections without boundary()==false
     - Id>=0 for all intersections without boundary()==true

     The way the Identifiers are attached to the grid may differ
     between the different grid implementations.

   */
  int boundaryId() const
  {
    return real_intersection_.boundaryId();
  }
#endif

  /** \brief index of the boundary segment within the macro grid
   *
   *  In many applications, special data needs to be attached to the boundary
   *  segments of the macro grid (e.g., a function selecting the boundary
   *  condition).
   *  Usually, this data is inherited by the children of the boundary segment.
   *
   *  In the DUNE framework, data is stored in arrays, addressed by an index,
   *  in this case the boundarySegmentIndex. The size of these arrays can be
   *  obtained by the Grid::numBoundarySegments.
   */
  size_t boundarySegmentIndex() const
  {
    return real_intersection_.boundarySegmentIndex();
  }

  //! @brief return true if intersection is shared with another element.
  bool neighbor() const
  {
    return true;
  }

  /*! @brief return EntityPointer to the Entity on the inside of this
     intersection. That is the Entity where we started this .
   */
  EntityPointer inside() const
  {
    return real_intersection_.inside();
  }

  /*! @brief return EntityPointer to the Entity on the outside of this
     intersection. That is the neighboring Entity.

     @warning Don't call this method if there is no neighboring Entity
     (neighbor() returns false). In this case the result is undefined.
   */
  EntityPointer outside() const
  {
    if (real_intersection_.neighbor())
      return real_intersection_.outside();
    else
      // do fancy stuff
      return real_intersection_.outside();
  }

  /*! @brief Return true if intersection is conforming.
   */
  bool conforming() const
  {
    return real_intersection_.conforming();
  }

  /** \brief geometrical information about this intersection in local
   *         coordinates of the inside() entity.
   *
   *  This method returns a Geometry object that provides a mapping from
   *  local coordinates of the intersection to local coordinates of the
   *  inside() entity.
   *
   *  \note Previously, the geometry was encapsulated in the intersection object
   *        and a const reference was returned.
   *
   *  \note The returned geometry object is guaranteed to remain valid until the
   *        grid is modified (or deleted).
   */
  LocalGeometry geometryInInside() const
  {
    return real_intersection_.geometryInInside();
  }

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
    if (real_intersection_.neighbor())
      return real_intersection_.geometryInOutside();
    else
      return real_intersection_.geometryInOutside(); // find periodic neighbor entity and corresponding intersection and
    // use geometryinInside for that Intersection
  }

  /** \brief geometrical information about the intersection in global coordinates.
   *
   *  This method returns a Geometry object that provides a mapping from
   *  local coordinates of the intersection to global (world) coordinates.
   *
   *  \note If the returned geometry has type <b>none</b> then only a limited set of features
   *        is availalbe for the geometry, i.e. center and volume.
   *
   *  \note Previously, the geometry was encapsulated in the intersection object
   *        and a const reference was returned.
   *
   *  \note The returned geometry object is guaranteed to remain valid until the
   *        grid is modified (or deleted).
   */
  Geometry geometry() const
  {
    return real_intersection_.geometry();
  }

  /** \brief obtain the type of reference element for this intersection */
  GeometryType type() const
  {
    return real_intersection_.type();
  }

  /** \brief Local index of codim 1 entity in the inside() entity where
   *         intersection is contained in
   *
   *  \note This method returns the face number with respect to the generic
   *        reference element.
   *
   *  \returns the index of the inside entity's face containing this
   *           intersection (with respect to the generic reference element)
   */
  int indexInInside() const
  {
    return real_intersection_.indexInInside();
  }

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
    if (real_intersection_.neighbor())
      return real_intersection_.geometryInOutside();
    else
      return 0; // find periodic neighbor entity and corresponding intersection and use indexInInside for that
    // Intersection
  }

  /*! @brief Return an outer normal (length not necessarily 1)

     The returned vector may depend on local position within the intersection.
   */
  GlobalCoordinate outerNormal(const LocalCoordinate& local) const
  {
    return real_intersection_.outerNormal(local);
  }

  /*! @brief return outer normal scaled with the integration element
        @copydoc Dune::Intersection::outerNormal
     The normal is scaled with the integration element of the intersection. This
        method is redundant but it may be more efficent to use this function
        rather than computing the integration element via geometry().
   */
  GlobalCoordinate integrationOuterNormal(const LocalCoordinate& local) const
  {
    return real_intersection_.integrationOuterNormal(local);
  }

  /*! @brief Return unit outer normal (length == 1)

     The returned vector may depend on the local position within the intersection.
     It is scaled to have unit length.
   */
  GlobalCoordinate unitOuterNormal(const LocalCoordinate& local) const
  {
    return real_intersection_.unitOuterNormal(local);
  }

  /*! @brief Return unit outer normal (length == 1)

     The returned vector is the normal at the center() of the
     intersection's geometry.
     It is scaled to have unit length.
   */
  GlobalCoordinate centerUnitOuterNormal() const
  {
    return real_intersection_.centerUnitOuterNormal();
  }

  //! Compares two intersections for equality.
  bool operator==(const PeriodicIntersection& other) const
  {
    return real_intersection_.equals(other.real_intersection_);
  }

  //! Compares two intersections for inequality.
  bool operator!=(const PeriodicIntersection& other) const
  {
    return !real_intersection_.equals(other.real_intersection_);
  }

protected:
  const RealIntersectionType& real_intersection_;
};

template <class GridViewImp>
class PeriodicIntersectionIterator
//     : IntersectionIterator < typename GridViewImp::Traits::Grid, PeriodicIntersectionIterator< GridViewImp >,
//     PeriodicIntersection< GridViewImp > >
{

  // type of underlying implementation, for internal use only
  typedef GridViewImp GridViewType;
  typedef PeriodicIntersectionIterator<GridViewType> ThisType;
  typedef PeriodicIntersection<GridViewType> PeriodicIntersectionType;
  typedef typename GridViewImp::IntersectionIterator RealIntersectionIteratorType;

public:
  /** \brief Type of Intersection this IntersectionIterator points to */
  typedef typename Dune::Intersection<typename GridViewType::Traits::Grid, PeriodicIntersectionType> Intersection;


  /** Copy Constructor from IntersectionIteratorImp */
  PeriodicIntersectionIterator(RealIntersectionIteratorType& real_intersection_iterator)
    : real_intersection_iterator_(real_intersection_iterator)
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
    return *(new Intersection(*(new PeriodicIntersectionType(*real_intersection_iterator_))));
  }

  /** \brief Pointer operator. */
  const Intersection* operator->() const
  {
    return (new Intersection(*(new PeriodicIntersectionType(*real_intersection_iterator_))));
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

  /** @brief Preincrement operator. Proceed to next intersection.*/
  ThisType& operator++()
  {
    ++real_intersection_iterator_;
    return *this;
  }

  /** @brief forward equality check to realIterator */
  bool equals(const PeriodicIntersectionIterator& rhs) const
  {
    return real_intersection_iterator_.equals(rhs.real_intersection_iterator_);
  }

protected:
  RealIntersectionIteratorType& real_intersection_iterator_;
};


template <class GridViewImp>
class PeriodicGridView;

template <class RealGridViewImp>
class PeriodicGridViewTraits
{
  typedef PeriodicGridView<RealGridViewImp> GridViewImp;
  typedef typename RealGridViewImp::Grid Grid;
  typedef typename RealGridViewImp::PartitionIterator PartitionIteratorType;

  /** \brief type of the index set */
  typedef typename Grid::Traits::LeafIndexSet IndexSet;

  /** \brief type of the intersection */
  typedef typename Grid::Traits::LeafIntersection Intersection;

  /** \brief type of the intersection iterator */
  typedef typename Grid::Traits::LeafIntersectionIterator IntersectionIterator;

  /** \brief type of the collective communication */
  typedef typename Grid::Traits::CollectiveCommunication CollectiveCommunication;

  template <int cd>
  struct Codim
  {
    typedef typename Grid::Traits::template Codim<cd>::template Partition<pitype>::LeafIterator Iterator;

    typedef typename Grid::Traits::template Codim<cd>::Entity Entity;
    typedef typename Grid::Traits::template Codim<cd>::EntityPointer EntityPointer;

    typedef typename Grid::template Codim<cd>::Geometry Geometry;
    typedef typename Grid::template Codim<cd>::LocalGeometry LocalGeometry;

    /** \brief Define types needed to iterate over entities of a given partition type */
    template <PartitionIteratorType pit>
    struct Partition
    {
      /** \brief iterator over a given codim and partition type */
      typedef typename Grid::template Codim<cd>::template Partition<pit>::LeafIterator Iterator;
    };
  };

  enum
  {
    conforming = Capabilities::isLeafwiseConforming<Grid>::v
  };
};
};


template <class GridViewImp>
class PeriodicGridView
//      : public GridView< typename GridViewImp::Traits >
{
  typedef PeriodicGridView<GridViewImp> ThisType;
  typedef GridView<typename GridViewImp::Traits> BaseType;
  typedef GridViewImp GridViewType;

public:
  typedef typename GridViewType::Grid Grid;
  typedef typename GridViewType::IndexSet IndexSet;
  typedef typename Dune::GeometryType GeometryType;
  typedef Dune::IntersectionIterator<Grid, PeriodicIntersectionIterator<GridViewType>,
                                     PeriodicIntersection<GridViewType>> IntersectionIterator;
  typedef typename GridViewType::CollectiveCommunication CollectiveCommunication;

  /** \brief constructor (engine concept) */
  PeriodicGridView(const GridViewType& real_grid_view)
      //      : BaseType(real_grid_view)
      : real_grid_view_(real_grid_view)
  {
  }

  /** \brief Copy constructor */
  PeriodicGridView(const ThisType& other)
      //      : BaseType( other.real_grid_view_ )
      : real_grid_view_(other.real_grid_view_)
  {
  }

  /** \brief assignment operator */
  ThisType& operator=(const ThisType& other)
  {
    real_grid_view_ = other.real_grid_view_;
    return *this;
  }

public:
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
  typename GridViewType::template Codim<cd>::Iterator begin() const
  {
    return real_grid_view_.template begin<cd>();
  }

  /** \brief obtain end iterator for this view */
  template <int cd>
  typename GridViewType::template Codim<cd>::Iterator end() const
  {
    return real_grid_view_.template end<cd>();
  }

  /** \brief obtain begin iterator for this view */
  template <int cd, PartitionIteratorType pitype>
  typename GridViewType::template Codim<cd>::template Partition<pitype>::Iterator begin() const
  {
    return real_grid_view_.template begin<cd, pitype>();
  }

  /** \brief obtain end iterator for this view */
  template <int cd, PartitionIteratorType pitype>
  typename GridViewType::template Codim<cd>::template Partition<pitype>::Iterator end() const
  {
    return real_grid_view_.template end<cd, pitype>();
  }

  /** \brief obtain begin intersection iterator with respect to this view */
  IntersectionIterator ibegin(const typename GridViewType::template Codim<0>::Entity& entity) const
  {
    IntersectionIterator periodic_intersection_iterator = *(new IntersectionIterator(real_grid_view_.ibegin(entity)));
    return periodic_intersection_iterator;
  }

  /** \brief obtain end intersection iterator with respect to this view */
  IntersectionIterator iend(const typename GridViewType::template Codim<0>::Entity& entity) const
  {
    IntersectionIterator periodic_intersection_iterator = *(new IntersectionIterator(real_grid_view_.ibegin(entity)));
    return periodic_intersection_iterator;
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

protected:
  const GridViewType& real_grid_view_;
};

} // namespace Grid
} // namespace Stuff
} // namespace Dune

#endif // #ifndef DUNE_STUFF_GRID_PERIODICVIEW_HH
