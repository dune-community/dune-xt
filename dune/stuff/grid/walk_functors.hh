// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_WALK_FUNCTORS_HH
#define DUNE_STUFF_WALK_FUNCTORS_HH

#include "walk.hh"

#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/grid/information.hh>

namespace Dune {
namespace Stuff {
namespace Grid {

//! GridWalk functor that refines all entitites above given volume
template <class GridType>
struct MaximumEntityVolumeRefineFunctor
{
  MaximumEntityVolumeRefineFunctor(GridType& grid, double volume, double factor)
    : threshold_volume_(volume * factor)
    , grid_(grid)
  {
  }

  template <class Entity>
  void operator()(const Entity& ent, const int /*ent_idx*/)
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
  MaximumEntityVolumeRefineFunctor<GridType> f(grid, unrefined_min_volume, size_factor);
  while (true) {
    grid.preAdapt();
    GridWalk<View>(view).walkCodim0(f);
    if (!grid.adapt())
      break;
    grid.postAdapt();
    std::cout << DSG::Dimensions<GridType>()(grid);
  }
} // EnforceMaximumEntityVolume

/** \brief Functor for a \ref GridWalk calculating minima and maxima of entities' coordinates
 **/
template <class EntityType>
struct MinMaxCoordinateFunctor
{
  typedef typename EntityType::Geometry EntityGeometryType;
  typedef typename EntityGeometryType::ctype ctype;
  typedef FieldVector<ctype, EntityGeometryType::coorddimension> VectorType;
  MinMaxCoordinateFunctor()
    : minima_(VectorType(std::numeric_limits<ctype>::max()))
    , maxima_(VectorType(std::numeric_limits<ctype>::min()))
  {
  }

  void operator()(const EntityType& ent, const int)
  {
    const typename EntityType::Geometry& geo = ent.geometry();
    for (auto i : DSC::valueRange(geo.corners())) {
      for (auto k : valueRange(EntityGeometryType::coorddimension)) {
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


#endif // DUNE_STUFF_WALK_FUNCTORS_HH
