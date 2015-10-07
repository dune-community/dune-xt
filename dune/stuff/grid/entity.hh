// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_GRID_ENTITY_HH
#define DUNE_STUFF_GRID_ENTITY_HH

#include <dune/geometry/referenceelements.hh>

#if HAVE_DUNE_GRID
#include <dune/grid/common/entity.hh>
#include <dune/grid/common/gridview.hh>
#endif

#include <dune/stuff/aliases.hh>
#include <dune/stuff/common/print.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/common/string.hh>
#include <dune/stuff/common/type_utils.hh>

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
  for (auto ii : DSC::valueRange(geometry.corners()))
    out << prefix + "  "
        << "corner " + Common::toString(ii) << " = " << geometry.corner(ii) << "\n";
} // ... printEntity(...)

#if HAVE_DUNE_GRID

template <int codim, int worlddim, class GridImp, template <int, int, class> class EntityImp>
double entity_diameter(const Dune::Entity<codim, worlddim, GridImp, EntityImp>& entity)
{
  auto max_dist        = std::numeric_limits<typename GridImp::ctype>::min();
  const auto& geometry = entity.geometry();
  for (auto i : DSC::valueRange(geometry.corners())) {
    const auto xi = geometry.corner(i);
    for (auto j : DSC::valueRange(i + 1, geometry.corners())) {
      auto xj = geometry.corner(j);
      xj -= xi;
      max_dist = std::max(max_dist, xj.two_norm());
    }
  }
  return max_dist;
} // entity_diameter

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

#endif // HAVE_DUNE_GRID

} // namespace Grid
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_GRID_ENTITY_HH
