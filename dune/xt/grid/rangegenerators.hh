// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2017 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)

#ifndef DUNE_XT_GRID_RANGEGENERATORS_HH
#define DUNE_XT_GRID_RANGEGENERATORS_HH

#include <dune/grid/common/rangegenerators.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/gridpart/common/gridpart.hh>
#endif

namespace Dune {
namespace XT {
namespace Grid {
namespace internal {


template <class GridLayerType>
class EntityRange
{
  const GridLayerType& layer_;

public:
  EntityRange(const GridLayerType& layer)
    : layer_(layer)
  {
  }

  auto begin() const -> decltype(layer_.template begin<0>())
  {
    return layer_.template begin<0>();
  }

  auto end() const -> decltype(layer_.template end<0>())
  {
    return layer_.template end<0>();
  }
}; // class EntityRange


template <class GridLayerType, class EntityType>
class IntersectionRange
{
  const GridLayerType& layer_;
  const EntityType& entity_;

public:
  IntersectionRange(const GridLayerType& layer, const EntityType& entity)
    : layer_(layer)
    , entity_(entity)
  {
  }

  auto begin() const -> decltype(layer_.ibegin(entity_))
  {
    return layer_.ibegin(entity_);
  }

  auto end() const -> decltype(layer_.iend(entity_))
  {
    return layer_.iend(entity_);
  }
}; // class IntersectionRange


} // namespace internal
} // namespace Grid
} // namespace XT

#if HAVE_DUNE_FEM


template <typename T>
inline auto elements(const Dune::Fem::GridPartInterface<T>& grid_part)
    -> decltype(XT::Grid::internal::EntityRange<Dune::Fem::GridPartInterface<T>>(grid_part))
{
  return XT::Grid::internal::EntityRange<Dune::Fem::GridPartInterface<T>>(grid_part);
}


template <typename T, typename Entity>
inline auto intersections(const Dune::Fem::GridPartInterface<T>& grid_part, const Entity& entity)
    -> decltype(XT::Grid::internal::IntersectionRange<Dune::Fem::GridPartInterface<T>, Entity>(grid_part, entity))
{
  return XT::Grid::internal::IntersectionRange<Dune::Fem::GridPartInterface<T>, Entity>(grid_part, entity);
}


#endif // HAVE_DUNE_FEM

} // namespace Dune

#endif // DUNE_XT_GRID_RANGEGENERATORS_HH
