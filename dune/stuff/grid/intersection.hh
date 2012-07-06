#ifndef DUNE_HELPER_TOOLS_GRID_INTERSECTION_HH
#define DUNE_FEMTOOLS_GRID_INTERSECTION_HH

// dune-common includes
#include <dune/common/fvector.hh>
#include <dune/common/static_assert.hh>

// dune-helper-tools includes
#include <dune/helper-tools/common/string.hh>

namespace Dune {

namespace HelperTools {

namespace Grid {

namespace Intersection {

/**
  \brief      prints some basic information about a Dune::Intersection, namely the number of its corners and the
              coordinates of those corners.
  \tparam     IntersectionType
              Dune::Intersection compatible
  \tparam     StreamType
              std::ostream compatible
  \param[in]  intersection
              Dune::Intersection, whose information should be printed
  \param[out] stream
              std::ostream, into which the information is printed
  **/
template <class IntersectionType, class StreamType = std::ostream>
void print(const IntersectionType& intersection, StreamType& stream = std::cout)
{
  typedef typename IntersectionType::Geometry GeometryType;

  typedef typename GeometryType::GlobalCoordinate GlobalPointType;

  const GeometryType& geometry = intersection.geometry();

  const int numCorners = geometry.corners();

  std::string prefix = "Dune::Intersection (" + Dune::HelperTools::Common::String::toString(numCorners) + " corner";
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
  const std::string whitespace = Dune::HelperTools::Common::String::whitespaceify(prefix);

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

template <class IntersectionType, class FieldType, int size>
bool contains(const IntersectionType& intersection, const Dune::FieldVector<FieldType, size>& globalPoint)
{
  dune_static_assert(size < 3,
                     "Dune::FemTools::Grid::Intersection::contains() not implemented for more than 2 dimension!");
  return false;
}

template <class IntersectionType, class FieldType>
bool contains(const IntersectionType& intersection, const Dune::FieldVector<FieldType, 1>& globalPoint)
{
  typedef Dune::FieldVector<FieldType, 1> GlobalPointType;

  typedef typename IntersectionType::Geometry GeometryType;

  const GeometryType& geometry = intersection.geometry();

  // get the only corner
  const GlobalPointType corner = geometry.corner(0);
  // check if the point is the corner
  if (corner == globalPoint) {
    return true;
  }

  return false;
} // end function contains

template <class IntersectionType, class FieldType>
bool contains(const IntersectionType& intersection, const Dune::FieldVector<FieldType, 2>& globalPoint)
{
  typedef Dune::FieldVector<FieldType, 2> GlobalPointType;

  typedef typename IntersectionType::Geometry GeometryType;

  const GeometryType& geometry = intersection.geometry();

  // get the two corners
  const GlobalPointType firstCorner  = geometry.corner(0);
  const GlobalPointType secondCorner = geometry.corner(1);

  // check, that point is on the line between the two points
  const FieldType x1 = (globalPoint[0] - firstCorner[0]) / (secondCorner[0] - firstCorner[0]);
  const FieldType x2 = (globalPoint[1] - firstCorner[1]) / (secondCorner[1] - firstCorner[1]);
  if (!(x1 > x2) && !(x1 < x2) && !(x1 < 0.0) && !(x1 > 1.0)) {
    return true;
  }

  return false;
} // end function contains

} // end namespace Intersection

} // end namespace Grid

} // end namespace HelperTools

} // end namespace Dune

#endif // DUNE_FEMTOOLS_GRID_INTERSECTION_HH
