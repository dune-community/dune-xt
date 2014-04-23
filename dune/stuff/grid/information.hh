// This file is part of the dune-stuff project:
//   https://users.dune-project.org/projects/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_GRID_INFORMATION_HH
#define DUNE_STUFF_GRID_INFORMATION_HH

#include <ostream>
#include <boost/format.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <dune/stuff/common/math.hh>
#include <dune/stuff/grid/intersection.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/grid/walk.hh>
#include <dune/stuff/aliases.hh>

namespace Dune {
namespace Stuff {
namespace Grid {

using namespace Dune::Stuff::Common;

struct Statistics
{
  int numberOfEntities;
  int numberOfIntersections;
  int numberOfInnerIntersections;
  int numberOfBoundaryIntersections;
  double maxGridWidth;
  template <class GridViewType>
  Statistics(const GridViewType& gridView)
    : numberOfEntities(gridView.size(0))
    , numberOfIntersections(0)
    , numberOfInnerIntersections(0)
    , numberOfBoundaryIntersections(0)
    , maxGridWidth(0)
  {
    for (const auto& entity : viewRange(gridView)) {
      for (const auto& intIt : intersectionRange(gridView, entity)) {
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
unsigned int maxNumberOfNeighbors(const GridViewType& gridView)
{
  unsigned int maxNeighbours = 0;
  for (const auto& entity : viewRange(gridView)) {
    unsigned int neighbours = 0;
    for (const auto& i : intersectionRange(gridView, entity)) {
      ++neighbours;
    }
    maxNeighbours = std::max(maxNeighbours, neighbours);
  }
  return maxNeighbours;
} // unsigned int maxNumberOfNeighbors(const GridPartType& gridPart)

//! Provide min/max coordinates for all space dimensions of a Grid (in the leafView)
template <class GridType>
struct Dimensions
{
  //! automatic running min/max
  typedef Dune::Stuff::Common::MinMaxAvg<typename GridType::ctype> MinMaxAvgType;
  typedef std::array<MinMaxAvgType, GridType::dimensionworld> CoordLimitsType;
  typedef typename GridType::template Codim<0>::Entity EntityType;
  CoordLimitsType coord_limits;
  MinMaxAvgType entity_volume;

  //! gridwalk functor that does the actual work for \ref GridDimensions
  class GridDimensionsFunctor
  {
    CoordLimitsType& coord_limits_;
    MinMaxAvgType& entity_volume_;

  public:
    GridDimensionsFunctor(CoordLimitsType& c, MinMaxAvgType& e)
      : coord_limits_(c)
      , entity_volume_(e)
    {
    }

    template <class Entity>
    void operator()(const Entity& ent, const int /*ent_idx*/)
    {
      const auto& geo = ent.geometry();
      entity_volume_(geo.volume());
      for (int i = 0; i < geo.corners(); ++i) {
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

  Dimensions(const GridType& grid)
  {
    typedef typename GridType::LeafGridView View;
    const auto& view = grid.leafView();
    GridDimensionsFunctor f(coord_limits, entity_volume);
    GridWalk<View>(view).walkCodim0(f);
  }

  Dimensions(const EntityType& entity)
  {
    GridDimensionsFunctor f(coord_limits, entity_volume);
    f(entity, 0);
  }
};

template <class GridType>
Dimensions<GridType> dimensions(const GridType& grid)
{
  return Dimensions<GridType>(grid);
}

template <class GridType>
Dimensions<GridType> dimensions(const typename GridType::template Codim<0>::Entity& entity)
{
  return Dimensions<GridType>(entity);
}

} // namespace Grid
} // end of namespace Stuff
} // namespace Dune

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
} // <<
#endif // DUNE_STUFF_GRID_INFORMATION_HH
