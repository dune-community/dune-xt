// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// The copyright lies with the authors of this file (see below).
// License: Dual licensed as  BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014, 2016)
//   Rene Milk       (2014 - 2015)
//   Tobias Leibner  (2015)

#ifndef DUNE_XT_GRID_WALKER_APPLY_ON_HH
#define DUNE_XT_GRID_WALKER_APPLY_ON_HH

#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/entity.hh>
#include <dune/xt/grid/intersection.hh>

namespace Dune {
namespace XT {
namespace Grid {
namespace ApplyOn {

/**
 *  \brief Interface for functors to tell on which entity to apply.
 *
 *  The derived class has to provide a method with the following signature:
 *  \code
virtual bool apply_on(const GridViewType& grid_view, const EntityType& entity) const
{
  ...
}
\endcode
 */
template <class GridViewImp>
class WhichEntity
{
public:
  typedef GridViewImp GridViewType;
  typedef typename XT::Grid::Entity<GridViewType>::Type EntityType;

  virtual ~WhichEntity()
  {
  }

  virtual bool apply_on(const GridViewType& /*grid_view*/, const EntityType& /*entity*/) const = 0;
}; // class WhichEntity

/**
 *  \brief Selects all entities.
 */
template <class GridViewImp>
class AllEntities : public WhichEntity<GridViewImp>
{
  typedef WhichEntity<GridViewImp> BaseType;

public:
  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::EntityType EntityType;

  virtual bool apply_on(const GridViewType& /*grid_view*/, const EntityType& /*entity*/) const override final
  {
    return true;
  }
}; // class AllEntities

/**
 *  \brief Selects entities which have a boundary intersection.
 */
template <class GridViewImp>
class BoundaryEntities : public WhichEntity<GridViewImp>
{
  typedef WhichEntity<GridViewImp> BaseType;

public:
  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::EntityType EntityType;

  virtual bool apply_on(const GridViewType& /*grid_view*/, const EntityType& entity) const override final
  {
    return entity.hasBoundaryIntersections();
  }
}; // class BoundaryEntities

/**
 *  \brief Interface for functors to tell on which intersection to apply.
 *
 *  The derived class has to provide a method with the following signature:
 *  \code
virtual bool apply_on(const GridViewType& grid_view, const IntersectionType& intersection) const
{
  ...
}
\endcode
 */
template <class GridViewImp>
class WhichIntersection
{
public:
  typedef GridViewImp GridViewType;
  typedef typename XT::Grid::Intersection<GridViewType>::Type IntersectionType;

  virtual ~WhichIntersection<GridViewImp>()
  {
  }

  virtual bool apply_on(const GridViewType& /*grid_view*/, const IntersectionType& /*intersection*/) const = 0;
}; // class WhichIntersection< GridViewImp >

/**
 *  \brief Selects all intersections.
 */
template <class GridViewImp>
class AllIntersections : public WhichIntersection<GridViewImp>
{
  typedef WhichIntersection<GridViewImp> BaseType;

public:
  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::IntersectionType IntersectionType;

  virtual bool apply_on(const GridViewType& /*grid_view*/,
                        const IntersectionType& /*intersection*/) const override final
  {
    return true;
  }
}; // class AllIntersections

/**
 *  \brief Selects each inner intersection.
 *
 *  To decide if this in an inner intersection,
\code
intersection.neighbor() && !intersection.boundary()
\endcode
 *  is used.
 */
template <class GridViewImp>
class InnerIntersections : public WhichIntersection<GridViewImp>
{
  typedef WhichIntersection<GridViewImp> BaseType;

public:
  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::IntersectionType IntersectionType;

  virtual bool apply_on(const GridViewType& /*grid_view*/, const IntersectionType& intersection) const override final
  {
    return intersection.neighbor() && !intersection.boundary();
  }
}; // class InnerIntersections

/**
 *  \brief Selects each inner intersection only once.
 *
 *  To decide if this in an inner intersection,
\code
intersection.neighbor() && !intersection.boundary()
\endcode
 *  is used, and true is returned, if the index of the inside() entity is smaller than the index of the outside()
 *  entity.
 */
template <class GridViewImp>
class InnerIntersectionsPrimally : public WhichIntersection<GridViewImp>
{
  typedef WhichIntersection<GridViewImp> BaseType;

public:
  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::IntersectionType IntersectionType;

  virtual bool apply_on(const GridViewType& grid_view, const IntersectionType& intersection) const override final
  {
    if (intersection.neighbor() && !intersection.boundary()) {
      const auto insideEntityPtr    = intersection.inside();
      const auto& insideEntity      = *insideEntityPtr;
      const auto outsideNeighborPtr = intersection.outside();
      const auto& outsideNeighbor = *outsideNeighborPtr;
      return grid_view.indexSet().index(insideEntity) < grid_view.indexSet().index(outsideNeighbor);
    } else
      return false;
  }
}; // class InnerIntersections

template <class GridViewImp>
class BoundaryIntersections : public WhichIntersection<GridViewImp>
{
  typedef WhichIntersection<GridViewImp> BaseType;

public:
  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::IntersectionType IntersectionType;

  virtual bool apply_on(const GridViewType& /*grid_view*/, const IntersectionType& intersection) const override final
  {
    return intersection.boundary();
  }
}; // class BoundaryIntersections

template <class GridViewImp>
class NonPeriodicBoundaryIntersections : public WhichIntersection<GridViewImp>
{
  typedef WhichIntersection<GridViewImp> BaseType;

public:
  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::IntersectionType IntersectionType;

  virtual bool apply_on(const GridViewType& /*grid_view*/, const IntersectionType& intersection) const override final
  {
    return intersection.boundary() && !intersection.neighbor();
  }
}; // class BoundaryIntersections

/**
 *  \brief Selects each periodic intersection.
 *
 *  To decide if this in an periodic intersection,
\code
intersection.neighbor() && intersection.boundary()
\endcode
 *  is used.
 */
template <class GridViewImp>
class PeriodicIntersections : public WhichIntersection<GridViewImp>
{
  typedef WhichIntersection<GridViewImp> BaseType;

public:
  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::IntersectionType IntersectionType;

  virtual bool apply_on(const GridViewType& /*grid_view*/, const IntersectionType& intersection) const override final
  {
    return intersection.neighbor() && intersection.boundary();
  }
}; // class PeriodicIntersections

template <class GridViewImp>
class FilteredIntersections : public WhichIntersection<GridViewImp>
{
  typedef WhichIntersection<GridViewImp> BaseType;

public:
  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::IntersectionType IntersectionType;
  typedef std::function<bool(const GridViewType&, const IntersectionType&)> FilterType;

  FilteredIntersections(FilterType filter)
    : filter_(filter)
  {
  }

  virtual bool apply_on(const GridViewType& grid_view, const IntersectionType& intersection) const override final
  {
    return filter_(grid_view, intersection);
  }

private:
  const FilterType filter_;
}; // class BoundaryIntersections

template <class GridViewImp>
class DirichletIntersections : public WhichIntersection<GridViewImp>
{
  typedef WhichIntersection<GridViewImp> BaseType;

public:
  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::IntersectionType IntersectionType;

  explicit DirichletIntersections(const BoundaryInfo<IntersectionType>& boundary_info)
    : boundary_info_(boundary_info)
  {
  }

  virtual bool apply_on(const GridViewType& /*grid_view*/, const IntersectionType& intersection) const override final
  {
    return boundary_info_.type(intersection) == DirichletBoundary();
  }

private:
  const BoundaryInfo<IntersectionType>& boundary_info_;
}; // class DirichletIntersections

template <class GridViewImp>
class NeumannIntersections : public WhichIntersection<GridViewImp>
{
  typedef WhichIntersection<GridViewImp> BaseType;

public:
  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::IntersectionType IntersectionType;

  explicit NeumannIntersections(const BoundaryInfo<IntersectionType>& boundary_info)
    : boundary_info_(boundary_info)
  {
  }

  virtual bool apply_on(const GridViewType& /*grid_view*/, const IntersectionType& intersection) const override final
  {
    return boundary_info_.type(intersection) == NeumannBoundary();
  }

private:
  const BoundaryInfo<IntersectionType>& boundary_info_;
}; // class NeumannIntersections

} // namespace ApplyOn
} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_WALKER_APPLY_ON_HH
