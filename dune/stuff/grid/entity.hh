#ifndef DUNE_STUFF_GRID_ENTITY_HH
#define DUNE_STUFF_GRID_ENTITY_HH

// dune-stuff includes
#include <dune/stuff/common/string.hh>

namespace Dune {
namespace Stuff {
namespace Grid {

template <class EntityType, class StreamType = std::ostream>
void printEntity(const EntityType& entity, StreamType& stream = std::cout, std::string prefix = "")
{
  const auto& geometry = entity.geometry();
  const int numCorners = geometry.corners();

  std::string header = "Dune::Entity (" + DSC::toString(numCorners) + " corner";
  if (numCorners != 1) {
    header += "s";
  }
  header += "): ";
  const unsigned int missing = 32 - header.size();
  if (missing > 0) {
    for (unsigned int i = 0; i < missing; ++i) {
      header += " ";
    }
  }
  const std::string whitespace = DSC::whitespaceify(header);

  stream << prefix << header << "[ (";
  for (int i = 0; i < numCorners; ++i) {
    const auto& corner = geometry.corner(i);
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

template <class GridImp, template <int, int, class> class EntityImp>
double geometryDiameter(const Dune::Entity<0, 2, GridImp, EntityImp>& entity)
{
  const auto end = entity.ileafend();
  double factor = 1.0;
  for (auto it = entity.ileafbegin(); it != end; ++it) {
    const auto& intersection = *it;
    factor *= intersection.geometry().volume();
  }
  return factor / (2.0 * entity.geometry().volume());
} // geometryDiameter

template <class GridImp, template <int, int, class> class EntityImp>
double geometryDiameter(const Dune::Entity<0, 3, GridImp, EntityImp>& entity)
{
  DUNE_THROW(Dune::NotImplementedError, "");
} // geometryDiameter

} // namespace Grid
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_GRID_ENTITY_HH

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
