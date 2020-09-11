// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012, 2014 - 2018)
//   René Fritze     (2012 - 2019)
//   Sven Kaulmann   (2014)
//   Tobias Leibner  (2014, 2016 - 2018, 2020)

#ifndef DUNE_XT_GRID_INFORMATION_HH
#define DUNE_XT_GRID_INFORMATION_HH

#include <ostream>

#include <boost/format.hpp>
#include <boost/range/adaptor/reversed.hpp>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/common/math.hh>
#include <dune/xt/common/ranges.hh>
#include <dune/xt/grid/entity.hh>
#include <dune/xt/grid/intersection.hh>
#include <dune/xt/grid/walker.hh>
#include <dune/xt/grid/functors/interfaces.hh>
#include <dune/xt/grid/type_traits.hh>

namespace Dune {
namespace XT {
namespace Grid {


struct Statistics
{
  size_t numberOfEntities;
  size_t numberOfIntersections;
  size_t numberOfInnerIntersections;
  size_t numberOfBoundaryIntersections;
  double maxGridWidth;
  template <class GridLayerType>
  Statistics(const GridLayerType& grid_layer)
    : numberOfEntities(grid_layer.size(0))
    , numberOfIntersections(0)
    , numberOfInnerIntersections(0)
    , numberOfBoundaryIntersections(0)
    , maxGridWidth(0)
  {
    for (auto&& entity : elements(grid_layer)) {
      for (auto&& intersection : intersections(grid_layer, entity)) {
        ++numberOfIntersections;
        maxGridWidth = std::max(intersection.geometry().volume(), maxGridWidth);
        // if we are inside the grid
        numberOfInnerIntersections += (intersection.neighbor() && !intersection.boundary());
        // if we are on the boundary of the grid
        numberOfBoundaryIntersections += (!intersection.neighbor() && intersection.boundary());
      }
    }
  }
};

/** \brief grid statistic output to given stream
 */
template <class GridLayerType>
void print_info(const GridLayerType& grid_layer, std::ostream& out)
{
  const Statistics st(grid_layer);
  out << "found " << st.numberOfEntities << " entities," << std::endl;
  out << "found " << st.numberOfIntersections << " intersections," << std::endl;
  out << "      " << st.numberOfInnerIntersections << " intersections inside and" << std::endl;
  out << "      " << st.numberOfBoundaryIntersections << " intersections on the boundary." << std::endl;
  out << "      maxGridWidth is " << st.maxGridWidth << std::endl;
} // ... print_info(...)

/**
 * \attention Not optimal, does a whole grid walk!
 **/
template <class GridLayerType>
size_t max_number_of_neighbors(const GridLayerType& grid_layer)
{
  size_t maxNeighbours = 0;
  for (auto&& entity : elements(grid_layer)) {
    size_t neighbours = 0;
    for (auto&& intersection : intersections(grid_layer, entity)) {
      (void)intersection; // silence unused variable warning
      ++neighbours;
    }
    maxNeighbours = std::max(maxNeighbours, neighbours);
  }
  return maxNeighbours;
} // ... max_number_of_neighbors(...)

//! Provide min/max coordinates for all space dimensions of a grid layer
template <class GridLayerType>
struct Dimensions
{
  static_assert(is_layer<GridLayerType>::value, "");
  typedef extract_grid_t<GridLayerType> GridType;
  //! automatic running min/max
  typedef Dune::XT::Common::MinMaxAvg<typename GridType::ctype> MinMaxAvgType;
  typedef std::array<MinMaxAvgType, GridType::dimensionworld> CoordLimitsType;
  typedef typename GridType::template Codim<0>::Entity ElementType;
  using BoundingBoxType = std::array<FieldVector<typename GridType::ctype, GridType::dimensionworld>, 2>;
  CoordLimitsType coord_limits;
  MinMaxAvgType entity_volume;
  MinMaxAvgType entity_width;

  //! gridwalk functor that does the actual work for \ref GridDimensions
  class GridDimensionsFunctor : public ElementFunctor<GridLayerType>
  {
    CoordLimitsType& coord_limits_;
    MinMaxAvgType& entity_volume_;
    MinMaxAvgType& entity_width_;

  public:
    GridDimensionsFunctor(CoordLimitsType& c, MinMaxAvgType& e, MinMaxAvgType& w)
      : coord_limits_(c)
      , entity_volume_(e)
      , entity_width_(w)
    {}

    virtual void apply_local(const ElementType& element) override
    {
      const auto& geo = element.geometry();
      entity_volume_(geo.volume());
      entity_width_(entity_diameter(element));
      for (auto i : Common::value_range(geo.corners())) {
        const auto& corner(geo.corner(i));
        for (size_t k = 0; k < GridType::dimensionworld; ++k)
          coord_limits_[k](corner[k]);
      }
    } // ()

    ElementFunctor<GridLayerType>* copy() override
    {
      return new GridDimensionsFunctor(*this);
    }
  };

  double volume_relation() const
  {
    return entity_volume.min() != 0.0 ? entity_volume.max() / entity_volume.min() : -1;
  }

  Dimensions(const GridLayerType& grid_layer)
  {
    GridDimensionsFunctor f(coord_limits, entity_volume, entity_width);
    Walker<GridLayerType> gw(grid_layer);
    gw.append(f);
    gw.walk();
  }

  Dimensions(const ElementType& entity)
  {
    GridDimensionsFunctor f(coord_limits, entity_volume, entity_width);
    f.apply_local(entity);
  }

  XT::Common::FieldVector<typename GridType::ctype, GridType::dimensionworld> view_center() const
  {
    XT::Common::FieldVector<typename GridType::ctype, GridType::dimensionworld> center;
    size_t idx = 0;
    for (auto&& axis_min_max : coord_limits) {
      center[idx++] = axis_min_max.average();
    }
    return center;
  }

  BoundingBoxType bounding_box() const
  {
    BoundingBoxType box;
    for (auto i : Common::value_range(GridType::dimensionworld)) {
      box[0][i] = coord_limits[i].min();
      box[1][i] = coord_limits[i].max();
    }
    return box;
  }
};

template <class GridType>
Dimensions<typename GridType::LeafGridViewType> dimensions(const GridType& grid)
{
  return Dimensions<typename GridType::LeafGridViewType>(grid.leafGridView());
}

template <class GridLayerType>
Dimensions<GridLayerType> dimensions(const GridLayerType& grid_layer)
{
  return Dimensions<GridLayerType>(grid_layer);
}

template <class GridLayerType>
Dimensions<GridLayerType> dimensions(const extract_entity_t<GridLayerType>& entity)
{
  return Dimensions<GridLayerType>(entity);
}

//! returns size() - overlap - ghosts
template <class GridType>
int parallel_size(const GridType& grid, int level, int codim)
{
  int size = grid.size(level, codim) - grid.overlapSize(level, codim) - grid.ghostSize(level, codim);
  return grid.comm().sum(size);
}

//! returns size() - overlap - ghosts
template <class GridType>
int parallel_size(const GridType& grid, int codim)
{
  return parallel_size(grid, grid.maxLevel(), codim);
}
} // namespace Grid
} // namespace XT
} // namespace Dune


template <class T>
inline std::ostream& operator<<(std::ostream& s, const Dune::XT::Grid::Dimensions<T>& d)
{
  for (size_t k = 0; k < T::dimensionworld; ++k) {
    const auto& mma = d.coord_limits[k];
    s << boost::format("x%d\tmin: %e\tavg: %e\tmax: %e\n") % k % mma.min() % mma.average() % mma.max();
  }
  s << boost::format("Entity vol min: %e\tavg: %e\tmax: %e\tQout: %e") % d.entity_volume.min()
           % d.entity_volume.average() % d.entity_volume.max() % d.volume_relation();
  s << std::endl;
  return s;
}


#endif // DUNE_XT_GRID_INFORMATION_HH
