// This file is part of the dune-stuff project:
//   https://users.dune-project.org/projects/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_GRID_ENTITY_HH
#define DUNE_STUFF_GRID_ENTITY_HH

#if HAVE_DUNE_GRID
#include <dune/grid/common/entity.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/gridview.hh>
#endif

#include <dune/stuff/common/string.hh>
#include <dune/stuff/common/print.hh>
#include <dune/stuff/aliases.hh>

namespace Dune {
namespace Stuff {
namespace Grid {

#if HAVE_DUNE_GRID


template <class GridPartOrViewType>
class Entity
{
  template <class GridViewType, bool is_view>
  struct Choose
  {
    typedef typename GridViewType::template Codim<0>::Entity Type;
  };

  template <class GridPartType>
  struct Choose<GridPartType, false>
  {
    typedef typename GridPartType::template Codim<0>::EntityType Type;
  };

  static const bool this_is_a_grid_view =
      std::is_base_of<GridView<typename GridPartOrViewType::Traits>, GridPartOrViewType>::value;

public:
  typedef typename Choose<GridPartOrViewType, this_is_a_grid_view>::Type Type;
}; // class Entity


#endif // HAVE_DUNE_GRID


template <class EntityType>
void printEntity(const EntityType& entity, const std::string name = Common::Typename<EntityType>::value(),
                 std::ostream& out = std::cout, const std::string prefix = "")
{
  if (!name.empty())
    out << prefix << name << ":\n";
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
  auto min_dist        = std::numeric_limits<typename GridImp::ctype>::max();
  const auto& geometry = entity.geometry();
  for (int i = 0; i < geometry.corners(); ++i) {
    const auto xi = geometry.corner(i);
    for (int j = i + 1; j < geometry.corners(); ++j) {
      auto xj = geometry.corner(j);
      xj -= xi;
      min_dist = std::min(min_dist, xj.two_norm());
    }
  }
  return min_dist;
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

template <int mydim, int cdim, class GridImp, template <int, int, class> class GeometryImp>
auto reference_element(const Dune::Geometry<mydim, cdim, GridImp, GeometryImp>& geometry)
    -> decltype(REFERENCE_ELEMENTS<typename GridImp::ctype, mydim>::general(geometry.type()))
{
  return REFERENCE_ELEMENTS<typename GridImp::ctype, mydim>::general(geometry.type());
}

#undef REFERENCE_ELEMENTS

} // namespace Grid
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_GRID_ENTITY_HH
