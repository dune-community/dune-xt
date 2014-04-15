#ifndef DUNE_STUFF_WALK_FUNCTORS_HH
#define DUNE_STUFF_WALK_FUNCTORS_HH

#include "walk.hh"

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
    for (int i = 0; i < geo.corners(); ++i) {
      for (int k = 0; k < int(EntityGeometryType::coorddimension); ++k) {
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

/** Copyright (c) 2012, Rene Milk, Felix Schindler
   * All rights reserved.
   *
   * Redistribution and use in source and binary forms, with or without
   * modification, are permitted provided that the following conditions are met:
   *
   * 1. Redistributions of source code must retain the above copyright notice, this
   *    list of conditions and the following disclaimer.
   * 2. Redistributions in binary form must reproduce the above copyright notice,
   *    this list of conditions and the following disclaimer in the documentation
   *    and/or other materials provided with the distribution.
   *
   * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
   * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
   * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
   * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
   * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
   * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
   * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
   * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
   *
   * The views and conclusions contained in the software and documentation are those
   * of the authors and should not be interpreted as representing official policies,
   * either expressed or implied, of the FreeBSD Project.
   **/
