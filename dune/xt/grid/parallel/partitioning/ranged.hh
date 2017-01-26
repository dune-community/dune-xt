#ifndef DUNE_XT_GRID_PARALLEL_PARTITIONING_RANGED_HH
#define DUNE_XT_GRID_PARALLEL_PARTITIONING_RANGED_HH

#include <dune/common/version.hh>
#if DUNE_VERSION_NEWER(DUNE_GRID, 3, 9) // EXADUNE
#include <dune/grid/utility/partitioning/ranged.hh>
#else
#include <dune/xt/grid/parallel/partitioning/ranged-internal.hh>
#endif

#endif // DUNE_XT_GRID_PARALLEL_PARTITIONING_RANGED_HH
