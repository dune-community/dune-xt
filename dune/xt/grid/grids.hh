// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// The copyright lies with the authors of this file (see below).
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016)

#ifndef DUNE_XT_GRID_GRIDS_HH
#define DUNE_XT_GRID_GRIDS_HH

#if HAVE_ALBERTA // clang-format off
# include <dune/xt/common/disable_warnings.hh>
#   include <dune/grid/albertagrid.hh>
# include <dune/xt/common/reenable_warnings.hh>
#endif
#include <dune/grid/yaspgrid.hh>

#if HAVE_DUNE_ALUGRID
# include <dune/alugrid/grid.hh>
#endif

#if HAVE_DUNE_SPGRID
# include <dune/grid/spgrid.hh>
#endif

#if HAVE_DUNE_UGGRID
# include <dune/grid/uggrid.hh>
#endif // clang-format on

#endif // DUNE_XT_GRID_GRIDS_HH
