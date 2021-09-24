// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   René Fritze     (2017 - 2019)
//   Tobias Leibner  (2017 - 2018, 2020)

#ifndef DUNE_XT_GRID_PARALLEL_PARTITIONING_RANGED_HH
#define DUNE_XT_GRID_PARALLEL_PARTITIONING_RANGED_HH

#include <dune/common/version.hh>
#if DUNE_VERSION_GTE(DUNE_GRID, 3, 9) // EXADUNE
#  include <dune/grid/utility/partitioning/ranged.hh>
#else
#  include <dune/xt/grid/parallel/partitioning/ranged-internal.hh>
namespace Dune {
template <class GridView, int codim>
using RangedPartitioning = Dune::XT::Grid::RangedPartitioning<GridView, codim, All_Partition>;
}
#endif

#endif // DUNE_XT_GRID_PARALLEL_PARTITIONING_RANGED_HH
