#ifndef DUNE_XT_GRID_GRIDS_HH
#define DUNE_XT_GRID_GRIDS_HH

#if HAVE_ALBERTA // clang-format off
# include <dune/xt/common/disable_warnings.hh>
#   include <dune/grid/albertagrid.hh>
# include <dune/xt/common/reenable_warnings.hh>
#endif
#if HAVE_ALUGRID
# include <dune/grid/alugrid.hh>
#endif
# include <dune/grid/sgrid.hh>
#include <dune/grid/yaspgrid.hh>

#if HAVE_DUNE_ALUGRID
# include <dune/alugrid/grid.hh>
#endif

#if HAVE_DUNE_SPGRID
# include <dune/grid/spgrid.hh>
#endif // clang-format on

#endif // DUNE_XT_GRID_GRIDS_HH
