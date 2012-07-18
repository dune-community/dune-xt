#ifndef DUNE_HELPER_TOOLS_GRID_ENTITY_HH
#define DUNE_HELPER_TOOLS_GRID_ENTITY_HH

// dune-stuff includes
#include <dune/stuff/common/string.hh>

namespace Dune {

namespace Stuff {

namespace Grid {

namespace Entity {

template <class EntityType>
void print(const EntityType& entity, std::ostream& stream = std::cout)
{
  typedef typename EntityType::Geometry GeometryType;

  typedef typename GeometryType::GlobalCoordinate GlobalPointType;

  const GeometryType& geometry = entity.geometry();

  const int numCorners = geometry.corners();

  std::string prefix = "Dune::Entity (" + Dune::HelperTools::Common::String::toString(numCorners) + " corner";
  if (numCorners != 1) {
    prefix += "s";
  }
  prefix += "): ";
  const unsigned int missing = 32 - prefix.size();
  if (missing > 0) {
    for (unsigned int i = 0; i < missing; ++i) {
      prefix += " ";
    }
  }
  const std::string whitespace = Dune::FemTools::String::whitespaceify(prefix);

  stream << prefix << "[ (";
  for (int i = 0; i < numCorners; ++i) {
    const GlobalPointType corner = geometry.corner(i);
    for (unsigned int j = 0; j < corner.size; ++j) {
      stream << corner[j];
      if (j < corner.size - 1) {
        stream << ", ";
      }
    }
    stream << ")";
    if (i < geometry.corners() - 1) {
      stream << "," << std::endl << whitespace << "  (";
    } else {
      stream << " ]" << std::endl;
    }
  }

} // end function print

template <class GridImp, template <int, int, class> class EntityImp>
double geometryDiameter(const Dune::Entity<0, 2, GridImp, EntityImp>& entity)
{
  typedef Dune::Entity<0, 2, GridImp, EntityImp> EntityType;
  typedef typename EntityType::LeafIntersectionIterator IntersectionIteratorType;
  IntersectionIteratorType end = entity.ileafend();
  double factor = 1.0;
  for (IntersectionIteratorType it = entity.ileafbegin(); it != end; ++it) {
    const typename IntersectionIteratorType::Intersection& intersection = *it;
    factor *= intersection.geometry().volume();
  }
  return factor / (2.0 * entity.geometry().volume());
} // geometryDiameter

template <class GridImp, template <int, int, class> class EntityImp>
double geometryDiameter(const Dune::Entity<0, 3, GridImp, EntityImp>& entity)
{
  DUNE_THROW(Dune::Exception, "copypasta from 2D");
  typedef Dune::Entity<0, 3, GridImp, EntityImp> EntityType;
  typedef typename EntityType::LeafIntersectionIterator IntersectionIteratorType;
  IntersectionIteratorType end = entity.ileafend();
  double factor = 1.0;
  for (IntersectionIteratorType it = entity.ileafbegin(); it != end; ++it) {
    const typename IntersectionIteratorType::Intersection& intersection = *it;
    factor *= intersection.geometry().volume();
  }
  return factor / (2.0 * entity.geometry().volume());
} // geometryDiameter

} // namespace Entity

} // namespace Grid

} // namespace Stuff

} // namespace Dune

#endif // DUNE_HELPER_TOOLS_GRID_ENTITY_HH
