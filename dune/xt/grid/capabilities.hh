// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   René Fritze     (2018)

#ifndef DUNE_XT_GRID_CAPABILITIES_HH
#define DUNE_XT_GRID_CAPABILITIES_HH

#include <dune/common/version.hh>

#include "grids.hh"

namespace Dune {
namespace XT {
namespace Grid {


template <class G, bool anything = true>
struct has_boundary_id
{
#if DUNE_VERSION_GTE(DUNE_GRID, 2, 6)
  static const constexpr bool value = false;
#else
#  if defined(DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS) && DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS
  static const constexpr bool value = true;
#  else
  static const constexpr bool value = false;
#  endif
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
