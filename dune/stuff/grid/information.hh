// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_GRID_INFORMATION_HH
#define DUNE_STUFF_GRID_INFORMATION_HH

#include <ostream>

#include <boost/format.hpp>
#include <boost/range/adaptor/reversed.hpp>

#include <dune/common/unused.hh>

#if HAVE_DUNE_GRID
#include <dune/grid/common/gridview.hh>
#endif

#include <dune/stuff/common/math.hh>
#include <dune/stuff/grid/intersection.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/grid/walker.hh>
#include <dune/stuff/aliases.hh>
#include <dune/stuff/grid/entity.hh>
#include <dune/stuff/grid/walker/functors.hh>

namespace Dune {
namespace Stuff {
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
    for (const auto& entity : DSC::entityRange(gridView)) {
      for (const auto& intIt : DSC::intersectionRange(gridView, entity)) {
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
void printInfo(const GridViewType& gridView, std::ostream& out)
{
  const Statistics st(gridView);
  out << "found " << st.numberOfEntities << " entities," << std::endl;
  out << "found " << st.numberOfIntersections << " intersections," << std::endl;
  out << "      " << st.numberOfInnerIntersections << " intersections inside and" << std::endl;
  out << "      " << st.numberOfBoundaryIntersections << " intersections on the boundary." << std::endl;
  out << "      maxGridWidth is " << st.maxGridWidth << std::endl;
} // printGridInformation

/**
* \attention Not optimal, does a whole grid walk!
**/
template <class GridViewType>
size_t maxNumberOfNeighbors(const GridViewType& gridView)
{
  size_t maxNeighbours = 0;
  for (const auto& entity : DSC::entityRange(gridView)) {
    size_t neighbours = 0;
    for (const auto& DSC_UNUSED(i) : DSC::intersectionRange(gridView, entity)) {
      ++neighbours;
    }
    maxNeighbours = std::max(maxNeighbours, neighbours);
  }
  return maxNeighbours;
} // size_t maxNumberOfNeighbors(const GridPartType& gridPart)

#if HAVE_DUNE_GRID

//! Provide min/max coordinates for all space dimensions of a GridView
template <class GridViewType>
struct Dimensions
{
  static_assert(std::is_base_of<GridView<typename GridViewType::Traits>, GridViewType>::value,
                "GridViewType is no GridView");
  typedef typename GridViewType::Grid GridType;
  //! automatic running min/max
  typedef Dune::Stuff::Common::MinMaxAvg<typename GridType::ctype> MinMaxAvgType;
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
      for (auto i : DSC::valueRange(geo.corners())) {
        const auto& corner(geo.corner(i));
        for (size_t k = 0; k < GridType::dimensionworld; ++k)
          coord_limits_[k](corner[k]);
      }
    } // ()
  };

  double volumeRelation() const
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

#endif // HAVE_DUNE_GRID

} // namespace Grid
} // end of namespace Stuff
} // namespace Dune

#if HAVE_DUNE_GRID

template <class T>
inline std::ostream& operator<<(std::ostream& s, const DSG::Dimensions<T>& d)
{
  for (size_t k = 0; k < T::dimensionworld; ++k) {
    const auto& mma = d.coord_limits[k];
    s << boost::format("x%d\tmin: %e\tavg: %e\tmax: %e\n") % k % mma.min() % mma.average() % mma.max();
  }
  s << boost::format("Entity vol min: %e\tavg: %e\tmax: %e\tQout: %e") % d.entity_volume.min()
           % d.entity_volume.average() % d.entity_volume.max() % d.volumeRelation();
  s << std::endl;
  return s;
}

#endif // HAVE_DUNE_GRID

#endif // DUNE_STUFF_GRID_INFORMATION_HH
