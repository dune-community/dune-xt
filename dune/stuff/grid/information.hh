#ifndef DUNE_STUFF_GRID_INFORMATION_HH
#define DUNE_STUFF_GRID_INFORMATION_HH

#include <ostream>
#include <boost/format.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <dune/stuff/common/math.hh>
#include <dune/stuff/grid/intersection.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/grid/walk.hh>

namespace Dune {
namespace Stuff {
namespace Grid {
namespace Information {

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
      for (const auto& intIt : Dune::Stuff::intersectionRange(gridView, entity)) {
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
void print(const GridViewType& gridView, std::ostream& out)
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
  typedef Dune::Stuff::Common::Math::MinMaxAvg<typename GridType::ctype> MinMaxAvgType;
  typedef Dune::array<MinMaxAvgType, GridType::dimensionworld> CoordLimitsType;
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
    Dune::Stuff::GridWalk<View>(view).walkCodim0(f);
  }
};


} // namespace Information
} // namespace Grid
} // end of namespace Stuff
} // namespace Dune

template <class T>
inline std::ostream& operator<<(std::ostream& s, const Dune::Stuff::Grid::Information::Dimensions<T>& d)
{
  for (size_t k = 0; k < T::dimensionworld; ++k) {
    const typename Dune::Stuff::Grid::Information::Dimensions<T>::MinMaxAvgType& mma = d.coord_limits[k];
    s << boost::format("x%d\tmin: %e\tavg: %e\tmax: %e\n") % k % mma.min() % mma.average() % mma.max();
  }
  s << boost::format("Entity vol min: %e\tavg: %e\tmax: %e\tQout: %e") % d.entity_volume.min()
           % d.entity_volume.average() % d.entity_volume.max() % d.volumeRelation();
  s << std::endl;
  return s;
} // <<
#endif // DUNE_STUFF_GRID_INFORMATION_HH
