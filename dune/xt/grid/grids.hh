// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2020)
//   René Fritze     (2017 - 2020)
//   Tim Keil        (2018)
//   Tobias Leibner  (2016, 2020)

#ifndef DUNE_XT_GRID_GRIDS_HH
#define DUNE_XT_GRID_GRIDS_HH

#if HAVE_ALBERTA
#  include <dune/xt/common/disable_warnings.hh>
#  include <dune/grid/albertagrid.hh>
#  include <dune/xt/common/reenable_warnings.hh>
#endif

#if HAVE_DUNE_ALUGRID
#  include <dune/alugrid/grid.hh>
#endif

#if HAVE_DUNE_SPGRID
#  include <dune/grid/spgrid.hh>
#  include <dune/grid/spgrid/dgfparser.hh>
#endif

#if HAVE_DUNE_UGGRID || HAVE_UG
#  include <dune/grid/uggrid.hh>
#endif

#include <dune/grid/onedgrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/xt/common/tuple.hh>

// this is used by other headers
using ONED_1D = Dune::OneDGrid;
using YASP_1D_EQUIDISTANT_OFFSET = Dune::YaspGrid<1, Dune::EquidistantOffsetCoordinates<double, 1>>;
using YASP_2D_EQUIDISTANT_OFFSET = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>;
using YASP_3D_EQUIDISTANT_OFFSET = Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<double, 3>>;
using YASP_4D_EQUIDISTANT_OFFSET = Dune::YaspGrid<4, Dune::EquidistantOffsetCoordinates<double, 4>>;
#if HAVE_DUNE_ALUGRID
using ALU_2D_SIMPLEX_CONFORMING = Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>;
using ALU_2D_SIMPLEX_NONCONFORMING = Dune::ALUGrid<2, 2, Dune::simplex, Dune::nonconforming>;
using ALU_2D_CUBE = Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>;
using ALU_3D_SIMPLEX_CONFORMING = Dune::ALUGrid<3, 3, Dune::simplex, Dune::conforming>;
using ALU_3D_SIMPLEX_NONCONFORMING = Dune::ALUGrid<3, 3, Dune::simplex, Dune::nonconforming>;
using ALU_3D_CUBE = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming>;
#endif
#if HAVE_DUNE_UGGRID || HAVE_UG
using UG_2D = Dune::UGGrid<2>;
using UG_3D = Dune::UGGrid<3>;
#endif
#if HAVE_ALBERTA
using ALBERTA_2D = Dune::AlbertaGrid<2, 2>;
using ALBERTA_3D = Dune::AlbertaGrid<3, 3>;
#endif


namespace Dune::XT::Grid {
namespace internal {


// To give better error messages, required below.
template <size_t d>
class ThereIsNoSimplexGridAvailableInDimension
{};


} // namespace internal


/**
 * \note Alberta grids are missing here on purpose, these cannot be handled automatically very well.
 */
using Available1dGridTypes = std::tuple<ONED_1D, YASP_1D_EQUIDISTANT_OFFSET>;

/**
 * \note Alberta grids are missing here on purpose, these cannot be handled automatically very well.
 */
using Available2dGridTypes = std::tuple<YASP_2D_EQUIDISTANT_OFFSET
#if HAVE_DUNE_ALUGRID
                                        ,
                                        ALU_2D_SIMPLEX_CONFORMING,
                                        ALU_2D_SIMPLEX_NONCONFORMING,
                                        ALU_2D_CUBE
#endif
#if HAVE_DUNE_UGGRID || HAVE_UG
                                        ,
                                        UG_2D
#endif
                                        >;

/**
 * \note Alberta grids are missing here on purpose, these cannot be handled automatically very well.
 */
using Available3dGridTypes = std::tuple<YASP_3D_EQUIDISTANT_OFFSET
#if HAVE_DUNE_ALUGRID
                                        ,
                                        ALU_3D_SIMPLEX_CONFORMING,
                                        ALU_3D_SIMPLEX_NONCONFORMING,
                                        ALU_3D_CUBE
#endif
#if HAVE_DUNE_UGGRID || HAVE_UG
                                        ,
                                        UG_3D
#endif
                                        >;

/**
 * \note Alberta grids are missing here on purpose, these cannot be handled automatically very well.
 */
using AvailableGridTypes = Common::tuple_cat_t<Available1dGridTypes, Available2dGridTypes, Available3dGridTypes>;

} // namespace Dune::XT::Grid


using SIMPLEXGRID_1D = ONED_1D;
using SIMPLEXGRID_2D =
#if HAVE_DUNE_ALUGRID
    ALU_2D_SIMPLEX_CONFORMING;
#elif HAVE_DUNE_UGGRID || HAVE_UG
    UG_2D;
#else
    Dune::XT::Grid::internal::ThereIsNoSimplexGridAvailableInDimension<2>;
#endif
using SIMPLEXGRID_3D =
#if HAVE_DUNE_ALUGRID
    ALU_3D_SIMPLEX_CONFORMING;
#elif HAVE_DUNE_UGGRID || HAVE_UG
    UG_3D;
#else
    Dune::XT::Grid::internal::ThereIsNoSimplexGridAvailableInDimension<3>;
#endif


using CUBEGRID_1D = ONED_1D;
using CUBEGRID_2D = YASP_2D_EQUIDISTANT_OFFSET;
using CUBEGRID_3D = YASP_3D_EQUIDISTANT_OFFSET;


#if HAVE_DUNE_ALUGRID || HAVE_DUNE_UGGRID || HAVE_UG
#  define SIMPLEXGRID_2D_AVAILABLE 1
#  define SIMPLEXGRID_3D_AVAILABLE 1
#else
#  define SIMPLEXGRID_2D_AVAILABLE 0
#  define SIMPLEXGRID_3D_AVAILABLE 0
#endif


using GRID_1D = ONED_1D;
using GRID_2D =
#if SIMPLEXGRID_2D_AVAILABLE
    SIMPLEXGRID_2D;
#else
    YASP_2D_EQUIDISTANT_OFFSET;
#endif
using GRID_3D =
#if SIMPLEXGRID_3D_AVAILABLE
    SIMPLEXGRID_3D;
#else
    YASP_3D_EQUIDISTANT_OFFSET;
#endif


#endif // DUNE_XT_GRID_GRIDS_HH
