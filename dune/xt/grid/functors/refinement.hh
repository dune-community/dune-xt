// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012, 2014 - 2018)
//   René Fritze     (2012, 2015 - 2016, 2018)

#ifndef DUNE_XT_GRID_FUNCTORS_REFINEMENT_HH
#define DUNE_XT_GRID_FUNCTORS_REFINEMENT_HH

#include <dune/xt/common/ranges.hh>
#include <dune/xt/grid/information.hh>

#include "interfaces.hh"

namespace Dune {
namespace XT {
namespace Grid {


//! GridWalk functor that refines all entitites above given volume
template <class GridViewType>
struct MaximumEntityVolumeRefineFunctor : public ElementFunctor<GridViewType>
{
  typedef ElementFunctor<GridViewType> BaseType;
  typedef typename GridViewType::GridType GridType;

  MaximumEntityVolumeRefineFunctor(GridType& grid, double volume, double factor)
    : threshold_volume_(volume * factor)
    , grid_(grid)
  {}

  void apply_local(const typename BaseType::ElementType& element) override final
  {
    const double volume = element.geometry().volume();
    if (volume > threshold_volume_)
      grid_.mark(1, element);
  }

  ElementFunctor<GridViewType>* copy() override
  {
    return new MaximumEntityVolumeRefineFunctor<GridViewType>(*this);
  }

  const double threshold_volume_;
  GridType& grid_;
};


//! refine entities until all have volume < size_factor * unrefined_minimum_volume
template <class GridType>
void enforce_maximum_entity_volume(GridType& grid, const double size_factor)
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
} // enforce_maximum_entity_volume


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_FUNCTORS_REFINEMENT_HH
