#ifndef DUNE_XT_GRID_CAPABILITIES_HH
#define DUNE_XT_GRID_CAPABILITIES_HH

#include "grids.hh"

namespace Dune {
namespace XT {
namespace Grid {


template <class G, bool anything = true>
struct has_boundary_id
{
#if defined(DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS) && DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS
  static const constexpr bool value = true;
#else
  static const constexpr bool value = false;
#endif
};


#if HAVE_DUNE_UGGRID || HAVE_UG

template <int dim, bool anything>
struct has_boundary_id<UGGrid<dim>, anything>
{
  static const constexpr bool value = false;
};

#endif // HAVE_DUNE_UGGRID || HAVE_UG


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_CAPABILITIES_HH
