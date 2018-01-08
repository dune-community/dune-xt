// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2017 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2017)
//   Tobias Leibner  (2016)

#ifndef DUNE_XT_GRID_GRIDS_HH
#define DUNE_XT_GRID_GRIDS_HH

#if HAVE_ALBERTA
#include <dune/xt/common/disable_warnings.hh>
#include <dune/grid/albertagrid.hh>
#include <dune/xt/common/reenable_warnings.hh>
#endif

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#if HAVE_DUNE_SPGRID
#include <dune/grid/spgrid.hh>
#include <dune/grid/spgrid/dgfparser.hh>
#endif

#if HAVE_DUNE_UGGRID || HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

#include <dune/grid/onedgrid.hh>
#include <dune/grid/yaspgrid.hh>


// this is used by other headers
typedef Dune::OneDGrid ONED_1D;
typedef Dune::YaspGrid<1, Dune::EquidistantOffsetCoordinates<double, 1>> YASP_1D_EQUIDISTANT_OFFSET;
typedef Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>> YASP_2D_EQUIDISTANT_OFFSET;
typedef Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<double, 3>> YASP_3D_EQUIDISTANT_OFFSET;
#if HAVE_DUNE_ALUGRID
typedef Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming> ALU_2D_SIMPLEX_CONFORMING;
typedef Dune::ALUGrid<2, 2, Dune::simplex, Dune::nonconforming> ALU_2D_SIMPLEX_NONCONFORMING;
typedef Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming> ALU_2D_CUBE;
typedef Dune::ALUGrid<3, 3, Dune::simplex, Dune::conforming> ALU_3D_SIMPLEX_CONFORMING;
typedef Dune::ALUGrid<3, 3, Dune::simplex, Dune::nonconforming> ALU_3D_SIMPLEX_NONCONFORMING;
typedef Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming> ALU_3D_CUBE;
#endif
#if HAVE_DUNE_UGGRID || HAVE_UG
typedef Dune::UGGrid<2> UG_2D;
typedef Dune::UGGrid<3> UG_3D;
#endif
#if HAVE_ALBERTA
typedef Dune::AlbertaGrid<2, 2> ALBERTA_2D;
#endif


#endif // DUNE_XT_GRID_GRIDS_HH
