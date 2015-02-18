// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_GRID_WALKER_FUNCTORS_HH
#define DUNE_STUFF_GRID_WALKER_FUNCTORS_HH

// nothing here will compile w/o grid present
#if HAVE_DUNE_GRID

#include <dune/stuff/grid/entity.hh>
#include <dune/stuff/grid/intersection.hh>
#include <dune/stuff/grid/boundaryinfo.hh>

namespace Dune {
namespace Stuff {
namespace Grid {
namespace Functor {


template <class GridViewImp>
class Codim0
{
public:
  typedef GridViewImp GridViewType;
  typedef typename Stuff::Grid::Entity<GridViewType>::Type EntityType;

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
  typedef typename Stuff::Grid::Entity<GridViewType>::Type EntityType;
  typedef typename Stuff::Grid::Intersection<GridViewType>::Type IntersectionType;

  virtual ~Codim1()
  {
  }

  virtual void prepare()
  {
  }

  virtual void apply_local(const IntersectionType& /*intersection*/, const EntityType& /*inside_entity*/,
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
  typedef typename Stuff::Grid::Entity<GridViewType>::Type EntityType;
  typedef typename Stuff::Grid::Intersection<GridViewType>::Type IntersectionType;

  virtual ~Codim0And1()
  {
  }

  virtual void prepare()
  {
  }

  virtual void apply_local(const EntityType& entity) = 0;

  virtual void apply_local(const IntersectionType& /*intersection*/, const EntityType& /*inside_entity*/,
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

  explicit DirichletDetector(const BoundaryInfoInterface<IntersectionType>& boundary_info)
    : boundary_info_(boundary_info)
    , found_(0)
  {
  }

  virtual ~DirichletDetector()
  {
  }

  virtual void apply_local(const IntersectionType& intersection, const EntityType& /*inside_entity*/,
                           const EntityType& /*outside_entity*/) override
  {
    if (boundary_info_.dirichlet(intersection))
      ++found_;
  }

  bool found() const
  {
    return found_ > 0;
  }

private:
  const BoundaryInfoInterface<IntersectionType>& boundary_info_;
  size_t found_;
}; // class DirichletDetector


} // namespace Functor
} // namespace Grid
} // namespace Stuff
} // namespace Dune

#endif // HAVE_DUNE_GRID

#endif // DUNE_STUFF_GRID_WALKER_FUNCTORS_HH
