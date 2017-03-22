// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2017 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014, 2016 - 2017)
//   Rene Milk       (2014 - 2016)
//   Tobias Leibner  (2015 - 2016)

#ifndef DUNE_XT_GRID_WALKER_APPLY_ON_HH
#define DUNE_XT_GRID_WALKER_APPLY_ON_HH

#include <dune/xt/common/memory.hh>

#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/entity.hh>
#include <dune/xt/grid/intersection.hh>

namespace Dune {
namespace XT {
namespace Grid {
namespace ApplyOn {


/**
 *  \brief Interface for functors to select on which entities to apply.
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

  virtual WhichEntity<GridViewImp>* copy() const = 0; // required for the python bindings

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
  using typename BaseType::GridViewType;
  typedef typename BaseType::EntityType EntityType;

  virtual WhichEntity<GridViewImp>* copy() const override final
  {
    return new AllEntities<GridViewImp>();
  }

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
  using typename BaseType::GridViewType;
  typedef typename BaseType::EntityType EntityType;

  virtual WhichEntity<GridViewImp>* copy() const override final
  {
    return new BoundaryEntities<GridViewImp>();
  }

  virtual bool apply_on(const GridViewType& /*grid_view*/, const EntityType& entity) const override final
  {
    return entity.hasBoundaryIntersections();
  }
}; // class BoundaryEntities


/**
 *  \brief Interface for functors to select on which intersections to apply.
 *
 *  \note  When deriving from this class, one can use the \sa internal::WhichIntersectionBase to provide a generic
 *         implementation of WhichIntersection::copy in two circumstances: (i) if one requires a BoundaryInfo, derive
 *         from internal::WhichIntersectionBase<..., true>, \sa DirichletIntersections; (ii) if the derived class does
 *         not have any members, derive from internal::WhichIntersectionBase<...>, \sa InnerIntersections. In all other
 *         circumstances, one has to manually implement copy, \sa FilteredIntersections.
 */
template <class GridViewImp>
class WhichIntersection
{
public:
  typedef GridViewImp GridViewType;
  using IntersectionType = extract_intersection_t<GridViewType>;

  virtual ~WhichIntersection<GridViewImp>()
  {
  }

  virtual WhichIntersection<GridViewImp>* copy() const = 0; // required for redirect lambdas, i.e., in python bindings

  virtual bool apply_on(const GridViewType& /*grid_view*/, const IntersectionType& /*intersection*/) const = 0;
};


namespace internal {


template <class GV, class Imp, bool ctor_with_boundary_info = false>
class WhichIntersectionBase : public WhichIntersection<GV>
{
public:
  virtual WhichIntersection<GV>* copy() const override final
  {
    return new Imp();
  }
};

template <class GV, class Imp>
class WhichIntersectionBase<GV, Imp, true> : public WhichIntersection<GV>
{
public:
  using typename WhichIntersection<GV>::IntersectionType;

  WhichIntersectionBase(const BoundaryInfo<IntersectionType>& boundary_info)
    : boundary_info_(boundary_info)
  {
  }

  virtual WhichIntersection<GV>* copy() const override final
  {
    return new Imp(boundary_info_.access());
  }

protected:
  Common::ConstStorageProvider<BoundaryInfo<IntersectionType>> boundary_info_;
};


} // namespace internal


/**
 *  \brief Selects all intersections.
 */
template <class GridViewImp>
class AllIntersections : public internal::WhichIntersectionBase<GridViewImp, AllIntersections<GridViewImp>>
{
  typedef WhichIntersection<GridViewImp> BaseType;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

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
class InnerIntersections : public internal::WhichIntersectionBase<GridViewImp, InnerIntersections<GridViewImp>>
{
  typedef WhichIntersection<GridViewImp> BaseType;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

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
class InnerIntersectionsPrimally
    : public internal::WhichIntersectionBase<GridViewImp, InnerIntersectionsPrimally<GridViewImp>>
{
  typedef WhichIntersection<GridViewImp> BaseType;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  virtual bool apply_on(const GridViewType& grid_view, const IntersectionType& intersection) const override final
  {
    if (intersection.neighbor() && !intersection.boundary()) {
      const auto insideEntity = intersection.inside();
      const auto outsideNeighbor = intersection.outside();
      return grid_view.indexSet().index(insideEntity) < grid_view.indexSet().index(outsideNeighbor);
    } else
      return false;
  }
}; // class InnerIntersections


template <class GridViewImp>
class BoundaryIntersections : public internal::WhichIntersectionBase<GridViewImp, BoundaryIntersections<GridViewImp>>
{
  typedef WhichIntersection<GridViewImp> BaseType;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  virtual bool apply_on(const GridViewType& /*grid_view*/, const IntersectionType& intersection) const override final
  {
    return intersection.boundary();
  }
}; // class BoundaryIntersections


template <class GridViewImp>
class NonPeriodicBoundaryIntersections
    : public internal::WhichIntersectionBase<GridViewImp, NonPeriodicBoundaryIntersections<GridViewImp>>
{
  typedef WhichIntersection<GridViewImp> BaseType;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

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
class PeriodicIntersections : public internal::WhichIntersectionBase<GridViewImp, PeriodicIntersections<GridViewImp>>
{
  typedef WhichIntersection<GridViewImp> BaseType;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  virtual bool apply_on(const GridViewType& /*grid_view*/, const IntersectionType& intersection) const override final
  {
    return intersection.neighbor() && intersection.boundary();
  }
}; // class PeriodicIntersections


template <class GridViewImp>
class PeriodicIntersectionsPrimally
    : public internal::WhichIntersectionBase<GridViewImp, PeriodicIntersectionsPrimally<GridViewImp>>
{
  typedef WhichIntersection<GridViewImp> BaseType;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  virtual bool apply_on(const GridViewType& grid_view, const IntersectionType& intersection) const override final
  {
    if (intersection.neighbor() && intersection.boundary()) {
      const auto insideEntity = intersection.inside();
      const auto outsideNeighbor = intersection.outside();
      return grid_view.indexSet().index(insideEntity) < grid_view.indexSet().index(outsideNeighbor);
    } else {
      return false;
    }
  }
}; // class PeriodicIntersectionsPrimally


template <class GridViewImp>
class FilteredIntersections : public WhichIntersection<GridViewImp>
{
  typedef WhichIntersection<GridViewImp> BaseType;
  typedef FilteredIntersections<GridViewImp> ThisType;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;
  typedef std::function<bool(const GridViewType&, const IntersectionType&)> FilterType;

  FilteredIntersections(FilterType filter)
    : filter_(filter)
  {
  }

  virtual WhichIntersection<GridViewImp>* copy() const override final
  {
    return new FilteredIntersections<GridViewImp>(filter_);
  }

  virtual bool apply_on(const GridViewType& grid_view, const IntersectionType& intersection) const override final
  {
    return filter_(grid_view, intersection);
  }

private:
  const FilterType filter_;
}; // class BoundaryIntersections


template <class GridViewImp>
class DirichletIntersections
    : public internal::WhichIntersectionBase<GridViewImp, DirichletIntersections<GridViewImp>, true>
{
  typedef internal::WhichIntersectionBase<GridViewImp, DirichletIntersections<GridViewImp>, true> BaseType;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  explicit DirichletIntersections(const BoundaryInfo<IntersectionType>& boundary_info)
    : BaseType(boundary_info)
  {
  }

  virtual bool apply_on(const GridViewType& /*grid_view*/, const IntersectionType& intersection) const override final
  {
    return boundary_info_.access().type(intersection) == DirichletBoundary();
  }

protected:
  using BaseType::boundary_info_;
}; // class DirichletIntersections


template <class GridViewImp>
class NeumannIntersections
    : public internal::WhichIntersectionBase<GridViewImp, NeumannIntersections<GridViewImp>, true>
{
  typedef internal::WhichIntersectionBase<GridViewImp, NeumannIntersections<GridViewImp>, true> BaseType;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  explicit NeumannIntersections(const BoundaryInfo<IntersectionType>& boundary_info)
    : BaseType(boundary_info)
  {
  }

  virtual bool apply_on(const GridViewType& /*grid_view*/, const IntersectionType& intersection) const override final
  {
    return boundary_info_.access().type(intersection) == NeumannBoundary();
  }

protected:
  using BaseType::boundary_info_;
}; // class NeumannIntersections


} // namespace ApplyOn
} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_WALKER_APPLY_ON_HH
