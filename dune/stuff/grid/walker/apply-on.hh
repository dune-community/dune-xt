// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_GRID_WALKER_APPLY_ON_HH
#define DUNE_STUFF_GRID_WALKER_APPLY_ON_HH

#include <dune/stuff/grid/entity.hh>
#include <dune/stuff/grid/intersection.hh>
#include <dune/stuff/grid/boundaryinfo.hh>

namespace Dune {
namespace Stuff {
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
  typedef typename Stuff::Grid::Entity<GridViewType>::Type EntityType;

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

  virtual bool apply_on(const GridViewType& /*grid_view*/, const EntityType& /*entity*/) const DS_OVERRIDE DS_FINAL
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

  virtual bool apply_on(const GridViewType& /*grid_view*/, const EntityType& entity) const DS_OVERRIDE DS_FINAL
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
  typedef typename Stuff::Grid::Intersection<GridViewType>::Type IntersectionType;

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
                        const IntersectionType& /*intersection*/) const DS_OVERRIDE DS_FINAL
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

  virtual bool apply_on(const GridViewType& /*grid_view*/,
                        const IntersectionType& intersection) const DS_OVERRIDE DS_FINAL
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

  virtual bool apply_on(const GridViewType& grid_view, const IntersectionType& intersection) const DS_OVERRIDE DS_FINAL
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

  virtual bool apply_on(const GridViewType& /*grid_view*/,
                        const IntersectionType& intersection) const DS_OVERRIDE DS_FINAL
  {
    return intersection.boundary();
  }
}; // class BoundaryIntersections


template <class GridViewImp>
class DirichletIntersections : public WhichIntersection<GridViewImp>
{
  typedef WhichIntersection<GridViewImp> BaseType;

public:
  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::IntersectionType IntersectionType;

  DirichletIntersections(const BoundaryInfoInterface<IntersectionType>& boundary_info)
    : boundary_info_(boundary_info)
  {
  }

  virtual bool apply_on(const GridViewType& /*grid_view*/,
                        const IntersectionType& intersection) const DS_OVERRIDE DS_FINAL
  {
    return boundary_info_.dirichlet(intersection);
  }

private:
  const BoundaryInfoInterface<IntersectionType>& boundary_info_;
}; // class DirichletIntersections


template <class GridViewImp>
class NeumannIntersections : public WhichIntersection<GridViewImp>
{
  typedef WhichIntersection<GridViewImp> BaseType;

public:
  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::IntersectionType IntersectionType;

  NeumannIntersections(const BoundaryInfoInterface<IntersectionType>& boundary_info)
    : boundary_info_(boundary_info)
  {
  }

  virtual bool apply_on(const GridViewType& /*grid_view*/,
                        const IntersectionType& intersection) const DS_OVERRIDE DS_FINAL
  {
    return boundary_info_.neumann(intersection);
  }

private:
  const BoundaryInfoInterface<IntersectionType>& boundary_info_;
}; // class NeumannIntersections


} // namespace ApplyOn
} // namespace Grid
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_GRID_WALKER_APPLY_ON_HH
