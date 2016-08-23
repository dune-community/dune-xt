// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// The copyright lies with the authors of this file (see below).
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014, 2016)
//   Rene Milk       (2014 - 2016)

#ifndef DUNE_XT_GRID_WALKER_FUNCTORS_HH
#define DUNE_XT_GRID_WALKER_FUNCTORS_HH

#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/entity.hh>
#include <dune/xt/grid/intersection.hh>

namespace Dune {
namespace XT {
namespace Grid {
namespace Functor {

template <class GridViewImp>
class Codim0
{
public:
  typedef GridViewImp GridViewType;
  typedef typename XT::Grid::Entity<GridViewType>::Type EntityType;

  virtual ~Codim0()
  {
  }

  virtual void prepare()
  {
  }

  virtual void apply_local(const EntityType& entity) = 0;

  virtual void finalize()
  {
  }
}; // class Codim0

template <class GridViewImp>
class Codim1
{
public:
  typedef GridViewImp GridViewType;
  typedef typename XT::Grid::Entity<GridViewType>::Type EntityType;
  typedef typename XT::Grid::Intersection<GridViewType>::Type IntersectionType;

  virtual ~Codim1()
  {
  }

  virtual void prepare()
  {
  }

  virtual void apply_local(const IntersectionType& /*intersection*/,
                           const EntityType& /*inside_entity*/,
                           const EntityType& /*outside_entity*/) = 0;

  virtual void finalize()
  {
  }
}; // class Codim1

template <class GridViewImp>
class Codim0And1
{
public:
  typedef GridViewImp GridViewType;
  typedef typename XT::Grid::Entity<GridViewType>::Type EntityType;
  typedef typename XT::Grid::Intersection<GridViewType>::Type IntersectionType;

  virtual ~Codim0And1()
  {
  }

  virtual void prepare()
  {
  }

  virtual void apply_local(const EntityType& entity) = 0;

  virtual void apply_local(const IntersectionType& /*intersection*/,
                           const EntityType& /*inside_entity*/,
                           const EntityType& /*outside_entity*/) = 0;

  virtual void finalize()
  {
  }
}; // class Codim0And1

template <class GridViewImp>
class DirichletDetector : public Codim1<GridViewImp>
{
  typedef Codim1<GridViewImp> BaseType;

public:
  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::IntersectionType IntersectionType;

  explicit DirichletDetector(const BoundaryInfo<IntersectionType>& boundary_info)
    : boundary_info_(boundary_info)
    , found_(0)
  {
  }

  virtual ~DirichletDetector()
  {
  }

  virtual void apply_local(const IntersectionType& intersection,
                           const EntityType& /*inside_entity*/,
                           const EntityType& /*outside_entity*/) override
  {
    if (boundary_info_.type(intersection) == DirichletBoundary())
      ++found_;
  }

  bool found() const
  {
    return found_ > 0;
  }

private:
  const BoundaryInfo<IntersectionType>& boundary_info_;
  size_t found_;
}; // class DirichletDetector

} // namespace Functor
} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_WALKER_FUNCTORS_HH
