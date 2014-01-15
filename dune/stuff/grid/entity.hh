#ifndef DUNE_STUFF_GRID_ENTITY_HH
#define DUNE_STUFF_GRID_ENTITY_HH

#if HAVE_DUNE_GRID
#include <dune/grid/common/entity.hh>
#include <dune/geometry/referenceelements.hh>
#endif

#include <dune/stuff/common/string.hh>
#include <dune/stuff/common/print.hh>
#include <dune/stuff/aliases.hh>

namespace Dune {
namespace Stuff {
namespace Grid {

template <class EntityType>
void printEntity(const EntityType& entity, std::ostream& out = std::cout, const std::string prefix = "")
{
  out << prefix << Common::Typename<EntityType>::value() << std::endl;
  const auto& geometry = entity.geometry();
  for (int ii = 0; ii < geometry.corners(); ++ii)
    Common::print(geometry.corner(ii), "corner " + Common::toString(ii), out, prefix + "  ");
} // ... printEntity(...)


#if HAVE_DUNE_GRID
template <class GridImp, template <int, int, class> class EntityImp>
double DUNE_DEPRECATED_MSG("use entityDiameter instead")
    geometryDiameter(const Dune::Entity<0, 2, GridImp, EntityImp>& entity)
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
double DUNE_DEPRECATED_MSG("use entityDiameter instead")
    geometryDiameter(const Dune::Entity<0, 3, GridImp, EntityImp>& /*entity*/)
{
  DUNE_THROW(Dune::NotImplemented, "geometryDiameter not implemented for dim 3");
} // geometryDiameter

template <int codim, int worlddim, class GridImp, template <int, int, class> class EntityImp>
double entityDiameter(const Dune::Entity<codim, worlddim, GridImp, EntityImp>& entity)
{
  const auto& geometry = entity.geometry();
  auto max_dist = std::numeric_limits<typename GridImp::ctype>::min();
  for (int i = 0; i < geometry.corners(); ++i) {
    const auto xi = geometry.corner(i);
    for (int j = i + 1; j < geometry.corners(); ++j) {
      auto xj = geometry.corner(j);
      xj -= xi;
      max_dist = std::max(max_dist, xj.two_norm());
    }
  }
  return max_dist;
} // geometryDiameter
#endif // HAVE_DUNE_GRID

#if DUNE_VERSION_NEWER(DUNE_GEOMETRY, 2, 3)
#define REFERENCE_ELEMENTS ReferenceElements
#else
#define REFERENCE_ELEMENTS GenericReferenceElements
#endif

template <int codim, int worlddim, class GridImp, template <int, int, class> class EntityImp>
auto reference_element(const Dune::Entity<codim, worlddim, GridImp, EntityImp>& entity)
    -> decltype(REFERENCE_ELEMENTS<typename GridImp::ctype, worlddim>::general(entity.geometry().type()))
{
  return REFERENCE_ELEMENTS<typename GridImp::ctype, worlddim>::general(entity.geometry().type());
}

#undef REFERENCE_ELEMENTS

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
