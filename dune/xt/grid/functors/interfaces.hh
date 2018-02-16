// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014, 2016 - 2017)
//   Rene Milk       (2014 - 2018)

#ifndef DUNE_XT_GRID_WALKER_FUNCTORS_HH
#define DUNE_XT_GRID_WALKER_FUNCTORS_HH

#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/entity.hh>
#include <dune/xt/grid/intersection.hh>

namespace Dune {
namespace XT {
namespace Grid {
namespace Functor {

template <class GridLayerImp>
class Codim0
{
public:
  typedef GridLayerImp GridLayerType;
  using EntityType = extract_entity_t<GridLayerType>;

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

template <class GridLayerImp, class ReturnImp>
class Codim0Return : public Codim0<GridLayerImp>
{
  typedef Codim0<GridLayerImp> BaseType;

public:
  typedef ReturnImp ReturnType;
  using typename BaseType::EntityType;

  virtual ~Codim0Return()
  {
  }

  virtual ReturnType compute_locally(const EntityType& entity) = 0;

  virtual ReturnType result() const = 0;
}; // class Codim0ReturnFunctor

template <class GridLayerImp>
class Codim1
{
public:
  typedef GridLayerImp GridLayerType;
  using EntityType = extract_entity_t<GridLayerType>;
  using IntersectionType = extract_intersection_t<GridLayerType>;

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

template <class GridLayerImp>
class Codim0And1
{
public:
  typedef GridLayerImp GridLayerType;
  using EntityType = extract_entity_t<GridLayerType>;
  using IntersectionType = extract_intersection_t<GridLayerType>;

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

template <class GridLayerImp>
class DirichletDetector : public Codim1<GridLayerImp>
{
  typedef Codim1<GridLayerImp> BaseType;

public:
  typedef typename BaseType::GridLayerType GridLayerType;
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
