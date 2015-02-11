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
class PeriodicIntersection : public GridViewImp::Intersection
{
  typedef PeriodicIntersection<GridViewImp> ThisType;
  typedef typename GridViewImp::Intersection BaseType;

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

template <class GridViewImp>
class PeriodicIntersectionIterator : public GridViewImp::IntersectionIterator
{

  // type of underlying implementation, for internal use only
  typedef GridViewImp GridViewType;
  typedef PeriodicIntersectionIterator<GridViewType> ThisType;
  typedef typename GridViewType::IntersectionIterator BaseType;
  typedef PeriodicIntersection<GridViewType> Intersection;

public:
  /** \brief Type of Intersection this IntersectionIterator points to */

  /** Copy Constructor from IntersectionIteratorImp */
  PeriodicIntersectionIterator(BaseType& real_intersection_iterator)
    : BaseType(real_intersection_iterator)
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
};

template <class GridViewImp>
class PeriodicGridView : public GridViewImp
{
  typedef PeriodicGridView<GridViewImp> ThisType;
  typedef GridViewImp BaseType;

public:
  typedef typename BaseType::Grid Grid;
  typedef typename BaseType::IndexSet IndexSet;
  typedef typename Dune::GeometryType GeometryType;
  typedef PeriodicIntersectionIterator<BaseType> IntersectionIterator;
  typedef PeriodicIntersection<BaseType> Intersection;
  typedef typename BaseType::CollectiveCommunication CollectiveCommunication;

  /** \brief constructor (engine concept) */
  PeriodicGridView(const BaseType& real_grid_view)
    : BaseType(real_grid_view)
  {
  }

  //  /** \brief Copy constructor */
  //  PeriodicGridView ( const ThisType &other )
  //    : BaseType( (BaseType&)other )
  //  {}

  /** \brief assignment operator */
  ThisType& operator=(const ThisType& other)
  {
    impl_ = other.impl_;
    return *this;
  }

public:
  /** \brief obtain begin intersection iterator with respect to this view */
  IntersectionIterator ibegin(const typename BaseType::template Codim<0>::Entity& entity) const
  {
    auto real_intersection_iterator      = ((BaseType*)this)->ibegin(entity);
    auto& real_intersection_iterator_ref = real_intersection_iterator;
    IntersectionIterator periodic_intersection_iterator(real_intersection_iterator_ref);
    return periodic_intersection_iterator;
  }

  /** \brief obtain end intersection iterator with respect to this view */
  IntersectionIterator iend(const typename BaseType::template Codim<0>::Entity& entity) const
  {
    auto real_intersection_iterator      = ((BaseType*)this)->iend(entity);
    auto& real_intersection_iterator_ref = real_intersection_iterator;
    IntersectionIterator periodic_intersection_iterator(real_intersection_iterator_ref);
    return periodic_intersection_iterator;
  }

protected:
  using BaseType::impl_;
};

} // namespace Grid
} // namespace Stuff
} // namespace Dune

#endif // #ifndef DUNE_STUFF_GRID_PERIODICVIEW_HH
