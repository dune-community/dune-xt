#ifndef DUNE_HELPER_TOOLS_GRID_ENTITY_HH
#define DUNE_HELPER_TOOLS_GRID_ENTITY_HH

// dune-helper-tools includes
#include <dune/helper-tools/common/string.hh>

namespace Dune {

namespace HelperTools {

namespace Grid {

namespace Entity {

template <class EntityType, class StreamType = std::ostream>
void print(const EntityType& entity, StreamType& stream = std::cout)
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

} // end namespace Entity

} // end namespace Grid

} // end namespace HelperTools

} // end namespace Dune

#endif // DUNE_HELPER_TOOLS_GRID_ENTITY_HH
