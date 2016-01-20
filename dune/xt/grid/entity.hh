// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2012 - 2016)
//   Rene Milk       (2012 - 2015)
//   Sven Kaulmann   (2014)
//   Tobias Leibner  (2014)

#ifndef DUNE_XT_GRID_ENTITY_HH
#define DUNE_XT_GRID_ENTITY_HH

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/entity.hh>
#include <dune/grid/common/gridview.hh>

#include <dune/xt/common/print.hh>
#include <dune/xt/common/ranges.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/common/type_utils.hh>

namespace Dune {
namespace XT {
namespace Grid {


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


template <class EntityType>
void printEntity(const EntityType& entity, const std::string name = Common::Typename<EntityType>::value(),
                 std::ostream& out = std::cout, const std::string prefix = "")
{
  if (!name.empty())
    out << prefix << name << ":\n";
  const auto& geometry = entity.geometry();
  for (auto ii : Common::value_range(geometry.corners()))
    out << prefix + "  "
        << "corner " + Common::to_string(ii) << " = " << geometry.corner(ii) << "\n";
} // ... printEntity(...)


template <int codim, int worlddim, class GridImp, template <int, int, class> class EntityImp>
double entity_diameter(const Dune::Entity<codim, worlddim, GridImp, EntityImp>& entity)
{
  auto max_dist        = std::numeric_limits<typename GridImp::ctype>::min();
  const auto& geometry = entity.geometry();
  for (auto i : Common::value_range(geometry.corners())) {
    const auto xi = geometry.corner(i);
    for (auto j : Common::value_range(i + 1, geometry.corners())) {
      auto xj = geometry.corner(j);
      xj -= xi;
      max_dist = std::max(max_dist, xj.two_norm());
    }
  }
  return max_dist;
} // entity_diameter


template <int codim, int worlddim, class GridImp, template <int, int, class> class EntityImp>
auto reference_element(const Dune::Entity<codim, worlddim, GridImp, EntityImp>& entity)
    -> decltype(ReferenceElements<typename GridImp::ctype, worlddim>::general(entity.geometry().type()))
{
  return ReferenceElements<typename GridImp::ctype, worlddim>::general(entity.geometry().type());
}

template <int mydim, int cdim, class GridImp, template <int, int, class> class GeometryImp>
auto reference_element(const Dune::Geometry<mydim, cdim, GridImp, GeometryImp>& geometry)
    -> decltype(ReferenceElements<typename GridImp::ctype, mydim>::general(geometry.type()))
{
  return ReferenceElements<typename GridImp::ctype, mydim>::general(geometry.type());
}


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_ENTITY_HH
