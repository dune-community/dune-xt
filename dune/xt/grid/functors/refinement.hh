// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012, 2014 - 2017)
//   Rene Milk       (2012, 2015 - 2016, 2018)

#ifndef DUNE_XT_WALK_FUNCTORS_HH
#define DUNE_XT_WALK_FUNCTORS_HH

#include <dune/xt/common/ranges.hh>

#include <dune/xt/grid/information.hh>

#include "walker.hh"
#include "walker/functors.hh"

namespace Dune {
namespace XT {
namespace Grid {

//! GridWalk functor that refines all entitites above given volume
template <class GridViewType>
struct MaximumEntityVolumeRefineFunctor : public Functor::Codim0<GridViewType>
{
  typedef Functor::Codim0<GridViewType> BaseType;
  typedef typename GridViewType::GridType GridType;
  MaximumEntityVolumeRefineFunctor(GridType& grid, double volume, double factor)
    : threshold_volume_(volume * factor)
    , grid_(grid)
  {
  }

  virtual void apply_local(const typename BaseType::EntityType& ent)
  {

    const double volume = ent.geometry().volume();
    if (volume > threshold_volume_)
      grid_.mark(1, ent);
  }

  const double threshold_volume_;
  GridType& grid_;
};

//! refine entities until all have volume < size_factor * unrefined_minimum_volume
template <class GridType>
void EnforceMaximumEntityVolume(GridType& grid, const double size_factor)
{
  using namespace Dune::XT;
  const typename Grid::Dimensions<GridType> unrefined_dimensions(grid);
  const double unrefined_min_volume = unrefined_dimensions.entity_volume.min();
  typedef typename GridType::LeafGridView View;
  View view = grid.leafView();
  MaximumEntityVolumeRefineFunctor<View> f(grid, unrefined_min_volume, size_factor);
  while (true) {
    grid.preAdapt();
    Walker<View> gw(view);
    gw.append(f);
    gw.walk();
    if (!grid.adapt())
      break;
    grid.postAdapt();
    std::cout << Grid::Dimensions<GridType>()(grid);
  }
} // EnforceMaximumEntityVolume

/** \brief Functor for a \ref GridWalk calculating minima and maxima of entities' coordinates
 **/
template <class GridViewType>
struct MinMaxCoordinateFunctor : public Functor::Codim0<GridViewType>
{
  typedef Functor::Codim0<GridViewType> BaseType;
  typedef typename BaseType::EntityType::Geometry EntityGeometryType;
  typedef typename EntityGeometryType::ctype ctype;
  typedef FieldVector<ctype, EntityGeometryType::coorddimension> VectorType;
  MinMaxCoordinateFunctor()
    : minima_(VectorType(std::numeric_limits<ctype>::max()))
    , maxima_(VectorType(std::numeric_limits<ctype>::min()))
  {
  }

  virtual void apply_local(const typename BaseType::EntityType& ent)
  {
    const auto& geo = ent.geometry();
    for (auto i : Common::value_range(geo.corners())) {
      for (auto k : Common::value_range(EntityGeometryType::coorddimension)) {
        minima_[k] = std::min(minima_[k], geo.corner(i)[k]);
        maxima_[k] = std::max(maxima_[k], geo.corner(i)[k]);
      }
    }
  }
  VectorType minima_;
  VectorType maxima_;
};

} // namespace Grid
} // namespace Stud
} // namespace Dune

#endif // DUNE_XT_WALK_FUNCTORS_HH
