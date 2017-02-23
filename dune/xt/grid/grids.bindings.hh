// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2017 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_XT_GRID_GRIDS_BINDINGS_HH
#define DUNE_XT_GRID_GRIDS_BINDINGS_HH

#include "grids.hh"

namespace Dune {


// this is used by other headers
typedef YaspGrid<2, EquidistantOffsetCoordinates<double, 2>> YASP_2D_EQUIDISTANT_OFFSET;
#if HAVE_ALUGRID || HAVE_DUNE_ALUGRID
typedef ALUGrid<2, 2, simplex, conforming> ALU_2D_SIMPLEX_CONFORMING;
#endif


} // namespace Dune

#endif // DUNE_XT_GRID_GRIDS_BINDINGS_HH
