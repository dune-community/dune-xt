#ifndef DUNE_STUFF_GRID_INFORMATION_HH
#define DUNE_STUFF_GRID_INFORMATION_HH

#include <ostream>
#include <boost/format.hpp>
#include <dune/stuff/common/math.hh>
#include <dune/stuff/grid/intersection.hh>
#include <dune/stuff/grid/walk.hh>

namespace Dune {
namespace Stuff {
namespace Grid {
namespace Information {

/** \brief grid statistic output to given stream
   * \todo not require a space to be passed
   */
template <class GridPartType, class DiscreteFunctionSpaceType>
void print(GridPartType& gridPart, DiscreteFunctionSpaceType& space, std::ostream& out)
{
  int numberOfEntities(0);
  int numberOfIntersections(0);
  int numberOfInnerIntersections(0);
  int numberOfBoundaryIntersections(0);
  double maxGridWidth(0.0);

  typedef typename GridPartType::GridType GridType;

  typedef typename GridType::template Codim<0>::Entity EntityType;

  typedef typename GridPartType::template Codim<0>::IteratorType EntityIteratorType;

  typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;

  EntityIteratorType entityItEndLog = space.end();
  for (EntityIteratorType entityItLog = space.begin(); entityItLog != entityItEndLog; ++entityItLog) {
    const EntityType& entity = *entityItLog;
    // count entities
    ++numberOfEntities;
    // walk the intersections
    IntersectionIteratorType intItEnd = gridPart.iend(entity);
    for (IntersectionIteratorType intIt = gridPart.ibegin(entity); intIt != intItEnd; ++intIt) {
      // count intersections
      ++numberOfIntersections;
      maxGridWidth = std::max(intIt->geometry().volume(), maxGridWidth);
      // if we are inside the grid
      if (intIt.neighbor() && !intIt.boundary()) {
        // count inner intersections
        ++numberOfInnerIntersections;
      }
      // if we are on the boundary of the grid
      if (!intIt.neighbor() && intIt.boundary()) {
        // count boundary intersections
        ++numberOfBoundaryIntersections;
      }
    }
  }
  out << "found " << numberOfEntities << " entities," << std::endl;
  out << "found " << numberOfIntersections << " intersections," << std::endl;
  out << "      " << numberOfInnerIntersections << " intersections inside and" << std::endl;
  out << "      " << numberOfBoundaryIntersections << " intersections on the boundary." << std::endl;
  out << "      maxGridWidth is " << maxGridWidth << std::endl;
} // printGridInformation

/**
* \attention Not optimal, does a whole grid walk!
**/
template <class GridPartType>
unsigned int maxNumberOfNeighbors(const GridPartType& gridPart)
{
  // some preparations
  unsigned int maxNeighbours = 0;
  unsigned int neighbours    = 0;
  typedef typename GridPartType::template Codim<0>::IteratorType EntityIteratorType;
  typedef typename GridPartType::template Codim<0>::EntityType EntityType;
  typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
  // walk over all entities
  const EntityIteratorType entityIteratorEnd = gridPart.template end<0>();
  for (EntityIteratorType entityIterator = gridPart.template begin<0>(); entityIterator != entityIteratorEnd;
       ++entityIterator) {
    const EntityType& entity = *entityIterator;
    neighbours               = 0;
    // walk over all neighbors
    const IntersectionIteratorType intersectionIteratorEnd = gridPart.iend(entity);
    for (IntersectionIteratorType intersectionIterator = gridPart.ibegin(entity);
         intersectionIterator != intersectionIteratorEnd;
         ++intersectionIterator) {
      ++neighbours;
    } // walk over all neighbors
    maxNeighbours = std::max(maxNeighbours, neighbours);
  } // walk over all entities
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
      typedef typename Entity::Geometry EntityGeometryType;
      typedef Dune::FieldVector<typename EntityGeometryType::ctype, EntityGeometryType::coorddimension> DomainType;
      const typename Entity::Geometry& geo = ent.geometry();
      entity_volume_(geo.volume());
      for (int i = 0; i < geo.corners(); ++i) {
        const DomainType& corner(geo.corner(i));
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
    const View& view = grid.leafView();
    GridDimensionsFunctor f(coord_limits, entity_volume);
    Dune::Stuff::Grid::Walk<View>(view).walkCodim0(f);
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
