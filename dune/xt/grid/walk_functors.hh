// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2012, 2014 - 2015)
//   Rene Milk       (2012, 2015)

#ifndef DUNE_XT_WALK_FUNCTORS_HH
#define DUNE_XT_WALK_FUNCTORS_HH

// nothing here will compile w/o grid present
#if HAVE_DUNE_GRID

#include "walker/functors.hh"
#include "walker.hh"

#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/grid/information.hh>

namespace Dune {
namespace Stuff {
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
  namespace DSG = Dune::Stuff::Grid;
  const typename DSG::Dimensions<GridType> unrefined_dimensions(grid);
  const double unrefined_min_volume = unrefined_dimensions.entity_volume.min();
  typedef typename GridType::LeafGridView View;
  View view = grid.leafView();
  MaximumEntityVolumeRefineFunctor<View> f(grid, unrefined_min_volume, size_factor);
  while (true) {
    grid.preAdapt();
    Walker<View> gw(view);
    gw.add(f);
    gw.walk();
    if (!grid.adapt())
      break;
    grid.postAdapt();
    std::cout << DSG::Dimensions<GridType>()(grid);
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
    for (auto i : DSC::valueRange(geo.corners())) {
      for (auto k : DSC::valueRange(EntityGeometryType::coorddimension)) {
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

#endif // HAVE_DUNE_GRID

#endif // DUNE_XT_WALK_FUNCTORS_HH
