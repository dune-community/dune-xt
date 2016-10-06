// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// The copyright lies with the authors of this file (see below).
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012, 2014 - 2016)
//   Rene Milk       (2012 - 2016)
//   Sven Kaulmann   (2014)
//   Tobias Leibner  (2014, 2016)

#ifndef DUNE_XT_GRID_INFORMATION_HH
#define DUNE_XT_GRID_INFORMATION_HH

#include <ostream>

#include <boost/format.hpp>
#include <boost/range/adaptor/reversed.hpp>

#include <dune/common/unused.hh>

#include <dune/grid/common/gridview.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/common/math.hh>
#include <dune/xt/common/ranges.hh>

#include <dune/xt/grid/entity.hh>
#include <dune/xt/grid/intersection.hh>
#include <dune/xt/grid/walker.hh>
#include <dune/xt/grid/walker/functors.hh>

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
  template <class GridViewType>
  Statistics(const GridViewType& gridView)
    : numberOfEntities(gridView.size(0))
    , numberOfIntersections(0)
    , numberOfInnerIntersections(0)
    , numberOfBoundaryIntersections(0)
    , maxGridWidth(0)
  {
    for (auto&& entity : elements(gridView)) {
      for (auto&& intIt : intersections(gridView, entity)) {
        ++numberOfIntersections;
        maxGridWidth = std::max(intIt.geometry().volume(), maxGridWidth);
        // if we are inside the grid
        numberOfInnerIntersections += (intIt.neighbor() && !intIt.boundary());
        // if we are on the boundary of the grid
        numberOfBoundaryIntersections += (!intIt.neighbor() && intIt.boundary());
      }
    }
  }
};

/** \brief grid statistic output to given stream
   */
template <class GridViewType>
void print_info(const GridViewType& gridView, std::ostream& out)
{
  const Statistics st(gridView);
  out << "found " << st.numberOfEntities << " entities," << std::endl;
  out << "found " << st.numberOfIntersections << " intersections," << std::endl;
  out << "      " << st.numberOfInnerIntersections << " intersections inside and" << std::endl;
  out << "      " << st.numberOfBoundaryIntersections << " intersections on the boundary." << std::endl;
  out << "      maxGridWidth is " << st.maxGridWidth << std::endl;
} // ... print_info(...)

/**
* \attention Not optimal, does a whole grid walk!
**/
template <class GridViewType>
size_t max_number_of_neighbors(const GridViewType& gridView)
{
  size_t maxNeighbours = 0;
  for (auto&& entity : elements(gridView)) {
    size_t neighbours = 0;
    for (auto&& DXTC_UNUSED(i) : intersections(gridView, entity)) {
      ++neighbours;
    }
    maxNeighbours = std::max(maxNeighbours, neighbours);
  }
  return maxNeighbours;
} // ... max_number_of_neighbors(...)

//! Provide min/max coordinates for all space dimensions of a GridView
template <class GridViewType>
struct Dimensions
{
  static_assert(std::is_base_of<GridView<typename GridViewType::Traits>, GridViewType>::value,
                "GridViewType is no GridView");
  typedef typename GridViewType::Grid GridType;
  //! automatic running min/max
  typedef Dune::XT::Common::MinMaxAvg<typename GridType::ctype> MinMaxAvgType;
  typedef std::array<MinMaxAvgType, GridType::dimensionworld> CoordLimitsType;
  typedef typename GridType::template Codim<0>::Entity EntityType;
  CoordLimitsType coord_limits;
  MinMaxAvgType entity_volume;
  MinMaxAvgType entity_width;

  //! gridwalk functor that does the actual work for \ref GridDimensions
  class GridDimensionsFunctor : public Functor::Codim0<GridViewType>
  {
    CoordLimitsType& coord_limits_;
    MinMaxAvgType& entity_volume_;
    MinMaxAvgType& entity_width_;

  public:
    GridDimensionsFunctor(CoordLimitsType& c, MinMaxAvgType& e, MinMaxAvgType& w)
      : coord_limits_(c)
      , entity_volume_(e)
      , entity_width_(w)
    {
    }

    virtual void apply_local(const EntityType& ent)
    {
      const auto& geo = ent.geometry();
      entity_volume_(geo.volume());
      entity_width_(entity_diameter(ent));
      for (auto i : Common::value_range(geo.corners())) {
        const auto& corner(geo.corner(i));
        for (size_t k = 0; k < GridType::dimensionworld; ++k)
          coord_limits_[k](corner[k]);
      }
    } // ()
  };

  double volume_relation() const
  {
    return entity_volume.min() != 0.0 ? entity_volume.max() / entity_volume.min() : -1;
  }

  Dimensions(const GridViewType& gridView)
  {
    GridDimensionsFunctor f(coord_limits, entity_volume, entity_width);
    Walker<GridViewType> gw(gridView);
    gw.add(f);
    gw.walk();
  }

  Dimensions(const EntityType& entity)
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
};

template <class GridType>
Dimensions<typename GridType::LeafGridViewType> dimensions(const GridType& grid)
{
  return Dimensions<typename GridType::LeafGridViewType>(grid.leafGridView());
}

template <class GridViewType>
Dimensions<GridViewType> dimensions(const GridViewType& gridView)
{
  return Dimensions<GridViewType>(gridView);
}

template <class GridViewType>
Dimensions<GridViewType> dimensions(const typename GridViewType::Grid::template Codim<0>::Entity& entity)
{
  return Dimensions<GridViewType>(entity);
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
