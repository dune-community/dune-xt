#ifndef DUNE_STUFF_GRID_INTERSECTION_HH
#define DUNE_STUFF_GRID_INTERSECTION_HH

#ifdef HAVE_CMAKE_CONFIG
#include "cmake_config.h"
#else
#include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#include <dune/common/fvector.hh>
#include <dune/common/static_assert.hh>

#include <dune/stuff/common/string.hh>
#include <dune/stuff/aliases.hh>

namespace Dune {
namespace Stuff {
namespace Grid {

/**
  \brief      prints some basic information about a Dune::Intersection, namely the number of its corners and the
              coordinates of those corners.
  \tparam     IntersectionType
              Dune::Intersection compatible
  \param[in]  intersection
              Dune::Intersection, whose information should be printed
  \param[out] stream
              std::ostream, into which the information is printed
  **/
template <class IntersectionType>
void printIntersection(const IntersectionType& intersection, std::ostream& stream = std::cout)
{
  const auto& geometry = intersection.geometry();
  const int numCorners = geometry.corners();
  std::string prefix = "Dune::Intersection (" + DSC::toString(numCorners) + " corner";
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
  const std::string whitespace = DSC::whitespaceify(prefix);

  stream << prefix << "[ (";
  for (int i = 0; i < numCorners; ++i) {
    const auto corner = geometry.corner(i);
    for (unsigned int j = 0; j < corner.size(); ++j) {
      stream << corner[j];
      if (j < corner.size() - 1) {
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
bool intersectionContains(const IntersectionType& /*intersection*/,
                          const Dune::FieldVector<FieldType, size>& /*globalPoint*/)
{
  dune_static_assert(size < 3,
                     "Dune::FemTools::Grid::Intersection::contains() not implemented for more than 2 dimension!");
  return false;
}

template <class IntersectionType, class FieldType>
bool intersectionContains(const IntersectionType& intersection, const Dune::FieldVector<FieldType, 1>& globalPoint)
{
  const auto& geometry = intersection.geometry();
  const auto corner    = geometry.corner(0);
  // check if the point is the corner
  return corner == globalPoint;
} // end function contains

template <class IntersectionType, class FieldType>
bool intersectionContains(const IntersectionType& intersection, const Dune::FieldVector<FieldType, 2>& globalPoint)
{
  const auto& geometry    = intersection.geometry();
  const auto firstCorner  = geometry.corner(0);
  const auto secondCorner = geometry.corner(1);

  // check, that point is on the line between the two points
  const FieldType x1 = (globalPoint[0] - firstCorner[0]) / (secondCorner[0] - firstCorner[0]);
  const FieldType x2 = (globalPoint[1] - firstCorner[1]) / (secondCorner[1] - firstCorner[1]);
  return (!(x1 > x2) && !(x1 < x2) && !(x1 < 0.0) && !(x1 > 1.0));
} // end function contains

} // end namespace Grid
} // end of namespace Stuff
} // end namespace Dune

#endif // DUNE_STUFF_GRID_INTERSECTION_HH
/** Copyright (c) 2012, Felix Albrecht
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
