#ifndef DUNE_STUFF_GRID_INTERSECTION_HH
#define DUNE_STUFF_GRID_INTERSECTION_HH

#ifdef HAVE_CMAKE_CONFIG
#include "cmake_config.h"
#else
#include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#include <dune/common/fvector.hh>
#include <dune/common/static_assert.hh>

#include <dune/geometry/referenceelements.hh>

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

/** Check whether a spatial point lies on an intersection.
*
* @param[in] intersection The intersection
* @param[in] globalPoint A Dune::FieldVector with the global coordinates of the point
* @return Returns true if the point lies on the intersection, false otherwise.
*/
template <class IntersectionType, class FieldType, int dim>
bool intersectionContains(const IntersectionType& intersection, const Dune::FieldVector<FieldType, dim>& globalPoint)
{
  // map global coordinates to local coordinates of the intersection
  const auto& intersectionGeometry = intersection.geometry();
  const auto& localPoint           = intersectionGeometry.local(globalPoint);

// get codim 1 reference element
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY, 2, 3)
  const auto& refElement = ReferenceElements<FieldType, dim - 1>::general(intersectionGeometry.type());
#else
  const auto& refElement = GenericReferenceElements<FieldType, dim - 1>::general(intersectionGeometry.type());
#endif
  // check whether reference element contains the local coordinates
  return refElement.checkInside(localPoint);
} // end function intersectionContains

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
